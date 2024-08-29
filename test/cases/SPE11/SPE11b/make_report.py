import numpy as np
import pathlib
import string
import re
import click
import pandas as pd
import matplotlib.pyplot as plt


def plot_df_over_time(df, col, filename, tscale=1):

    plt.figure()
    plt.plot(df["t [s]"] / tscale, df[col])
    plt.xscale("log")
    plt.savefig(filename)


YEAR = int(3.1536e7)  # seconds
HOUR = 3600  # seconds


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument("case", type=click.Choice(["a", "b", "c"], case_sensitive=False))
@click.option("--src", default="extracts")
@click.option("--out", default="report")
def cli(case, src, out):
    proceed(case, src, out)


def standard_files():
    dense_file_pattern = r"dense_data_([\d.e+-]+)\.csv"
    box_files = ["Box-A.csv", "Box-B.csv"]
    pop_files = ["POP1.csv", "POP2.csv"]
    convc_file = "convC.csv"
    seal_file = "SealTot.csv"
    # FixMe: the boundaryCO2_file is only for case b and c
    boundaryCO2_file = "BoundaryCO2.csv"
    return (
        dense_file_pattern,
        box_files,
        pop_files,
        convc_file,
        seal_file,
        boundaryCO2_file,
    )


def proceed(case, src, out, *files):
    src = pathlib.Path(src)
    out = pathlib.Path(out)

    if not files:
        files = standard_files()
    (
        dense_file_pattern,
        box_files,
        pop_files,
        convc_file,
        seal_file,
        boundaryCO2_file,
    ) = files
    pop_files = [src / f for f in pop_files]
    box_files = [src / f for f in box_files]
    convc_file = src / convc_file
    seal_file = src / seal_file
    boundaryCO2_file = src / boundaryCO2_file

    match case:
        case "a":
            unit1 = 1
            dt1 = 600
            out1 = out / "spe11a_time_series.csv"
            dt2 = range(0, 200 * HOUR, 1 * HOUR)
            dV2 = 0.01 * 0.01 * 1
            out2 = lambda s: out / f"spe11a_spatial_map_{s/HOUR:.0f}h.csv"
        case "b":
            unit1 = YEAR
            dt1 = 0.1
            out1 = out / "spe11b_time_series.csv"
            dt2 = [*range(1000 * YEAR, 2001 * YEAR, 5 * YEAR)]
            dV2 = 10 * 10 * 1
            out2 = lambda s: out / f"spe11b_spatial_map_{s/YEAR:.0f}y.csv"
        case "c":
            unit1 = YEAR
            dt1 = 0.1
            out1 = out / "spe11c_time_series.csv"
            dt2 = [
                *range(1000 * YEAR, 1050 * YEAR, 5 * YEAR),
                *range(1050 * YEAR, 1100 * YEAR, 25 * YEAR),
                *range(1100 * YEAR, 1500 * YEAR, 50 * YEAR),
                *range(1500 * YEAR, 2001 * YEAR, 100 * YEAR),
            ]
            dV2 = 50 * 50 * 10
            out2 = lambda s: out / f"spe11c_spatial_map_{s/YEAR:.0f}y.csv"
        case _:
            raise ValueError(f"case (={case!r}) must be 'a', 'b' or 'c'.")

    make_time_series(
        pop_files, box_files, convc_file, seal_file, boundaryCO2_file, dt1, unit1, out1
    )
    make_spatial_maps(src, dense_file_pattern, dt2, dV2, out2)
    make_performance_time_series()


def make_time_series(
    pop_files, box_files, convc_file, seal_file, boundaryCO2_file, dt, unit, out
):
    box_cols = lambda a: {
        f"mob{a}": f"mob{a} [kg]",
        f"imm{a}": f"imm{a} [kg]",
        f"diss{a}": f"diss{a} [kg]",
        f"seal{a}": f"seal{a} [kg]",
    }

    box = pd.concat(
        [
            pd.read_csv(file, index_col=False)
            .set_index("Time")[box_cols(a).keys()]
            .rename(columns=box_cols(a))
            for a, file in zip(string.ascii_uppercase, box_files)
        ],
        axis=1,
    )

    pop = pd.concat(
        [
            pd.read_csv(file, index_col=False)
            .set_index("Time")["node pressure"]
            .rename(f"p{i} [Pa]")
            for i, file in enumerate(pop_files, start=1)
        ],
        axis=1,
    )

    convC = pd.read_csv(convc_file, index_col=False).set_index("Time")["M_C [kg]"]

    seal = pd.read_csv(seal_file, index_col=False).set_index("Time")["SealTot [kg]"]

    boundaryCO2 = pd.read_csv(boundaryCO2_file, index_col=False).set_index("Time")[
        "BoundaryCO2 [kg]"
    ]

    df = pd.concat([pop, box, convC, seal, boundaryCO2], axis=1).reset_index(
        names="t [s]"
    )

    # eliminate first value (initial state)
    first_value = df["t [s]"].iloc[0]
    df = df[df["t [s]"] >= first_value + 1].reset_index(drop=True)
    # remove the first 1000 years (initialization phase)
    df["t [s]"] = df["t [s]"] - 1e3 * YEAR

    # time interpolation (get nearest value, could be more precise here)
    df_report = df.iloc[:0, :].copy()
    for t in range(0, int(1000.001 * YEAR), int(0.1 * YEAR)):
        df_line = (abs(df["t [s]"] - t)).idxmin()
        df_data = df.iloc[df_line].copy()
        df_data["t [s]"] = t
        df_report = df_report._append(df_data, ignore_index=True)

    df_report["BoundaryCO2 [kg]"] = 0.0

    # # strange, convert into years to reconvert into second
    # df["t [s]"] = df["t [s]"]/YEAR
    # first_value = df["t [s]"].iloc[0]
    # # remove the first 1000 years (initialization phase)
    # df["t [s]"] = df["t [s]"] - 1e3
    # #  convert into seconds and ???
    # fake_unit = "s"
    # one = pd.to_datetime(1, unit=fake_unit) - pd.to_datetime(1)
    # df["datetime"] = pd.to_datetime(df["t [s]"], unit=fake_unit)
    # # groupby = df.groupby(pd.Grouper(key="datetime", freq=dt * one))
    # # df = groupby.mean()  # FIXME devrait dépendre de l'espacement des temps
    # # df = df.dropna()
    # df["t [s]"] = (df.index - pd.to_datetime(1)) / one * unit
    # df = df.reset_index(drop=True)
    # # df = df[df["t [s]"] >= 0].reset_index(drop=True)
    # # first_value = df["t [s]"].iloc[0]
    # # df["t [s]"] = df["t [s]"] - first_value

    out.parent.mkdir(exist_ok=True)
    names = [" " + n for n in df.columns]
    names[0] = "#" + names[0]
    df_report.to_csv(out, index=False, float_format="% .6e", header=names)

    # plot_df_over_time(df_report, "p1 [Pa]", "pressure.png", tscale=YEAR)
    # plot_df_over_time(df_report, "mobA [kg]", "mobA.png", tscale=YEAR)


def make_spatial_maps(src, dense_file_pattern, dt, dV, out):
    fields = {
        "Points:0": "x [m]",
        "Points:2": "z [m]",
        "node pressure": "pressure [Pa]",
        "node gas saturation": "gas saturation [-]",
        "node CO2 mass fraction in liquid": "mass fraction of CO2 in liquid [-]",
        "node CO2 mass fraction in gas": "mass fraction of CO2 in vapor [-]",
        "node gas mass density": "phase mass density gas [kg/m3]",
        "node liquid mass density": "phase mass density water [kg/m3]",
        "node temperature": "temperature [C]",
    }

    load = lambda f: pd.read_csv(f)[fields.keys()].rename(columns=fields).astype(float)
    get_time = (
        # lambda f: float(m.group(1)) if (m := re.match(dense_file_pattern, f)) else None
        lambda f: float(f"{float(m.group(1))}")
        if (m := re.match(dense_file_pattern, f))
        else None
    )
    tags = {
        f: k
        for f in src.glob("dense_data_*.csv")
        if (k := get_time(f.name)) is not None
    }
    # print(f"Found {len(tags)} files matching the pattern.")

    # for file, time in tags.items():
    #     print(f"Processing file: {file} at time {time}")

    for time, df in iter_times(dt, tags.keys(), get_time, load):
        df["mass fraction of H2O in vapor [-]"] = (
            1 - df["mass fraction of CO2 in vapor [-]"]
        )
        df["total mass CO2 [kg]"] = dV * (
            df["gas saturation [-]"]
            * df["phase mass density gas [kg/m3]"]
            * df["mass fraction of CO2 in vapor [-]"]
            + (1 - df["gas saturation [-]"])
            * df["phase mass density water [kg/m3]"]
            * df["mass fraction of CO2 in liquid [-]"]
        )
        df = df[
            [
                "x [m]",
                "z [m]",
                "pressure [Pa]",
                "gas saturation [-]",
                "mass fraction of CO2 in liquid [-]",
                "mass fraction of H2O in vapor [-]",
                "phase mass density gas [kg/m3]",
                "phase mass density water [kg/m3]",
                "total mass CO2 [kg]",
                "temperature [C]",
            ]
        ]
        #
        names = [" " + n for n in df.columns]
        names[0] = "#" + names[0]
        df_filtred = df[(df["# x [m]"] != 0) & (df[" z [m]"] != 0)]
        df_filtred.to_csv(
            out(time - 1000 * YEAR), index=False, float_format="% .6e", header=names
        )


def iter_times(times, files, get_time, load, eps=1e16):
    """yield (time, df) for each time in times
    df is mixture of files close to time
    """
    # FIXME on fait du "nearest" mais on pourrait faire mieux
    for time in times:
        # print(f"Checking for time: {time/YEAR}")
        dist = lambda f: abs(get_time(f.name) - time)
        res = min(files, key=dist)
        if dist(res) <= eps:
            print(
                f"Selected file: {res} for time: {time/YEAR} with distance: {dist(res)/YEAR}"
            )
            # print(f"Time: {time/YEAR} - Using file: {res}")
            yield time, load(res)


# def iter_intervals(seq):
#     it = iter(seq)
#     start = next(it)
#     for stop in it:
#         yield start, stop
#         start = stop
# def iter_groups(seq, func, times):
#     for start, stop in iter_intervals(times):
#         group = [f for f in seq if start <= func(f) <= stop]
#         if group:
#             yield start, stop, group


def make_performance_time_series():
    ...  # TODO


if __name__ == "__main__":
    # proceed("b", "extracts", "report")
    cli()
