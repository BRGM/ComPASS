""" python pv_extract.py -h

How to use it with a isolated paraview installation:
    1. install paraview in a dedicated virtual environment:
        conda create -c conda-forge -n paraview-5.12 paraview=5.12
    2. run the script with `pvpython` from the env:
        conda run -n paraview-5.12 pvpython pv_extract.py "PATH/states.pvd" cells_rkt.vtu -o "RESULTS_PATH"

"""

import pathlib
import argparse


parser = argparse.ArgumentParser(
    description="Extracts sparse and dense data for SPE11. It uses paraview (and its special python interpreter `pvpython`) to load and interact with data.",
    usage="pvpython %(prog)s [-h] [-o DIR] SIMUL RKT",
    epilog=r"/!\ Only work if used with `pvpython` /!\.",
)
parser.add_argument("states", metavar="SIMUL", help="simulation file (states.pvd)")
parser.add_argument("rocktypes", metavar="RKT", help="rocktype file (cells_rkt.vtu)")
parser.add_argument(
    "-o", metavar="DIR", dest="output_dir", default="extracts", help="output directory"
)

args = parser.parse_args()
state_pvd_file = pathlib.Path(args.states).absolute()
rocktype_vtu_file = pathlib.Path(args.rocktypes).absolute()
extract_directory = pathlib.Path(args.output_dir).absolute()

print(f"simulation file (states.pvd)     : {state_pvd_file}")
print(f"rocktype file (cells_rkt.vtu)    : {rocktype_vtu_file}")
print(f"extraction directory (extracts/) : {extract_directory}")


#####################################################################
#####################################################################
#####################################################################

# ----------------------------------------------------------------
# paraview state begins here
# ----------------------------------------------------------------

# state file generated using paraview version 5.12.1
import paraview

paraview.compatibility.major = 5
paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
stateswithmasspvd = PVDReader(
    registrationName="states with mass.pvd", FileName=str(state_pvd_file)
)
stateswithmasspvd.CellArrays = []
stateswithmasspvd.PointArrays = [
    "node CO2 mass fraction in gas",
    "node CO2 mass fraction in liquid",
    "node CO2 solubility ratio",
    "node gas mass density",
    "node gas saturation",
    "node liquid mass density",
    "node pressure",
    "node temperature",
]

# *** WARNING ***:
# FOR SPE11B: since we have an extension in the Y direction we need to adapt the box geometries
# IF CHANGING THE MESH, CHANGE THE Y to coincide boxCgeometry.Origin = [X, Y, Z]
# create a new 'Plane'
boxCgeometry = Plane(registrationName="Box C geometry")
boxCgeometry.Origin = [3300.0, 2.0, 100.0]
boxCgeometry.Point1 = [7800.0, 2.0, 100.0]
boxCgeometry.Point2 = [3300.0, 2.0, 400.0]
boxCgeometry.XResolution = 900
boxCgeometry.YResolution = 60

# create a new 'Plane'
boxBgeometry = Plane(registrationName="Box B geometry")
boxBgeometry.Origin = [100.0, 2.0, 600.0]
boxBgeometry.Point1 = [3300.0, 2.0, 600.0]
boxBgeometry.Point2 = [100.0, 2.0, 1200.0]
boxBgeometry.XResolution = 640
boxBgeometry.YResolution = 120

# create a new 'Plane'
observationgeometry = Plane(registrationName="observation geometry")
observationgeometry.Origin = [0.0, 2.0, 0.0]
observationgeometry.Point1 = [8400.0, 2.0, 0.0]
observationgeometry.Point2 = [0.0, 2.0, 1200.0]
observationgeometry.XResolution = 840
observationgeometry.YResolution = 120

# create a new 'XML Unstructured Grid Reader'
rocktypes = XMLUnstructuredGridReader(
    registrationName="rocktypes", FileName=[str(rocktype_vtu_file)]
)

rocktypes.CellArrayStatus = ["rocktype"]
print(rocktypes.CellArrayStatus)
rocktypes.TimeArray = "None"

# create a new 'Append Attributes'
stateswithmassandrocktypes = AppendAttributes(
    registrationName="states with mass and rocktypes",
    Input=[stateswithmasspvd, rocktypes],
)

# create a new 'Resample With Dataset'
observation = ResampleWithDataset(
    registrationName="observation",
    SourceDataArrays=stateswithmassandrocktypes,
    DestinationMesh=observationgeometry,
)
observation.PassCellArrays = 1
observation.PassPointArrays = 1
observation.PassFieldArrays = 0
observation.ComputeTolerance = 0
observation.SnapToCellWithClosestPoint = 1
observation.CellLocator = "Static Cell Locator"

# create a new 'Resample With Dataset'
boxB = ResampleWithDataset(
    registrationName="Box B",
    SourceDataArrays=stateswithmassandrocktypes,
    DestinationMesh=boxBgeometry,
)
boxB.PassCellArrays = 1
boxB.PassPointArrays = 1
boxB.PassFieldArrays = 0
boxB.ComputeTolerance = 0
boxB.SnapToCellWithClosestPoint = 1
boxB.CellLocator = "Static Cell Locator"

# create a new 'Calculator'
massCO2ingas = Calculator(registrationName="mass CO2 in gas", Input=boxB)
massCO2ingas.ResultArrayName = "mCO2 gas"
massCO2ingas.Function = (
    '"node CO2 mass fraction in gas"*"node gas mass density"*"node gas saturation"'
)

# create a new 'Calculator'
massfreeCO2ingas = Calculator(
    registrationName="mass free CO2 in gas", Input=massCO2ingas
)
massfreeCO2ingas.ResultArrayName = "mobB"
massfreeCO2ingas.Function = '"mCO2 gas"*("node gas saturation">0.1)'

# create a new 'Calculator'
masstrappedCO2ingas = Calculator(
    registrationName="mass trapped CO2 in gas", Input=massfreeCO2ingas
)
masstrappedCO2ingas.ResultArrayName = "immB"
masstrappedCO2ingas.Function = '"mCO2 gas"*("node gas saturation"<=0.1)'

# create a new 'Calculator'
massCO2inliquid = Calculator(
    registrationName="mass CO2 in liquid", Input=masstrappedCO2ingas
)
massCO2inliquid.ResultArrayName = "dissB"
massCO2inliquid.Function = '"node CO2 mass fraction in liquid"*"node liquid mass density"*(1-"node gas saturation")'

# create a new 'Calculator'
massCO2insealrtk1 = Calculator(
    registrationName="mass CO2 in seal (rkt==1)", Input=[rocktypes, massCO2inliquid]
)
massCO2insealrtk1.ResultArrayName = "sealB"
# massCO2insealrtk1.Function = '("node CO2 mass fraction in liquid"*"node liquid mass density"*(1-"node gas saturation")+"node CO2 mass fraction in gas"*"node gas mass density"*"node gas saturation")*(rkt==1)'
massCO2insealrtk1.Function = '("mCO2 gas" + "dissB")*(rocktype==1)'

# create a new 'Integrate Variables'
integration = IntegrateVariables(
    registrationName="integration", Input=massCO2insealrtk1
)

# create a new 'Plot Data Over Time'
dataovertime = PlotDataOverTime(registrationName="data over time", Input=integration)
dataovertime.OnlyReportSelectionStatistics = 0

# create a new 'Probe Location'
pOP2 = ProbeLocation(
    registrationName="POP 2",
    Input=stateswithmassandrocktypes,
    ProbeType="Fixed Radius Point Source",
)

# init the 'Fixed Radius Point Source' selected for 'ProbeType'
pOP2.ProbeType.Center = [5100.0, 2.0, 1100.0]

# create a new 'Plot Data Over Time'
pressureovertime = PlotDataOverTime(registrationName="pressure over time", Input=pOP2)
pressureovertime.OnlyReportSelectionStatistics = 0

# create a new 'Plane'
boxAgeometry = Plane(registrationName="Box A geometry")
boxAgeometry.Origin = [3300.0, 2.0, 0.0]
boxAgeometry.Point1 = [8300.0, 2.0, 0.0]
boxAgeometry.Point2 = [3300.0, 2.0, 600.0]
boxAgeometry.XResolution = 1000
boxAgeometry.YResolution = 120

# create a new 'Resample With Dataset'
boxA = ResampleWithDataset(
    registrationName="Box A",
    SourceDataArrays=stateswithmassandrocktypes,
    DestinationMesh=boxAgeometry,
)
boxA.PassCellArrays = 1
boxA.PassPointArrays = 1
boxA.PassFieldArrays = 0
boxA.ComputeTolerance = 0
boxA.SnapToCellWithClosestPoint = 1
boxA.CellLocator = "Static Cell Locator"

# create a new 'Resample With Dataset'
boxC = ResampleWithDataset(
    registrationName="Box C",
    SourceDataArrays=stateswithmassandrocktypes,
    DestinationMesh=boxCgeometry,
)
boxC.PassCellArrays = 1
boxC.PassPointArrays = 1
boxC.PassFieldArrays = 0
boxC.ComputeTolerance = 0
boxC.SnapToCellWithClosestPoint = 1
boxC.CellLocator = "Static Cell Locator"

# create a new 'Probe Location'
pOP1 = ProbeLocation(
    registrationName="POP 1",
    Input=stateswithmassandrocktypes,
    ProbeType="Fixed Radius Point Source",
)

# init the 'Fixed Radius Point Source' selected for 'ProbeType'
pOP1.ProbeType.Center = [4500.0, 2.0, 500.0]

# create a new 'Plot Data Over Time'
pressureovertime_1 = PlotDataOverTime(registrationName="pressure over time", Input=pOP1)
pressureovertime_1.OnlyReportSelectionStatistics = 0

# create a new 'Calculator'
massC02ingas = Calculator(registrationName="mass C02 in gas", Input=boxA)
massC02ingas.ResultArrayName = "mCO2 gas"
massC02ingas.Function = (
    '"node CO2 mass fraction in gas"*"node gas mass density"*"node gas saturation"'
)

# create a new 'Calculator'
massfreeC02ingas = Calculator(
    registrationName="mass free C02 in gas", Input=massC02ingas
)
massfreeC02ingas.ResultArrayName = "mobA"
massfreeC02ingas.Function = '"mCO2 gas"*("node gas saturation">0.1)'

# create a new 'Calculator'
masstrappedCO2ingas_1 = Calculator(
    registrationName="mass trapped CO2 in gas", Input=massfreeC02ingas
)
masstrappedCO2ingas_1.ResultArrayName = "immA"
masstrappedCO2ingas_1.Function = '"mCO2 gas"*("node gas saturation"<=0.1)'

# create a new 'Calculator'
massCO2inliquid_1 = Calculator(
    registrationName="mass CO2 in liquid", Input=masstrappedCO2ingas_1
)
massCO2inliquid_1.ResultArrayName = "dissA"
massCO2inliquid_1.Function = '"node CO2 mass fraction in liquid"*"node liquid mass density"*(1-"node gas saturation")'

# create a new 'Calculator'
massCO2insealrtk1_1 = Calculator(
    registrationName="mass CO2 in seal (rkt==1)", Input=massCO2inliquid_1
)
massCO2insealrtk1_1.ResultArrayName = "sealA"
massCO2insealrtk1_1.Function = '("mCO2 gas" + "dissA")*(rocktype==1)'

# create a new 'Integrate Variables'
integration_1 = IntegrateVariables(
    registrationName="integration", Input=massCO2insealrtk1_1
)

# create a new 'Plot Data Over Time'
dataovertime_1 = PlotDataOverTime(
    registrationName="data over time", Input=integration_1
)
dataovertime_1.OnlyReportSelectionStatistics = 0

# create a new 'Cell Data to Point Data'
rktonnodes = CellDatatoPointData(
    registrationName="rkt on nodes", Input=stateswithmassandrocktypes
)
rktonnodes.ProcessAllArrays = 0
rktonnodes.CellDataArraytoprocess = ["rocktype"]

# create a new 'Calculator'
sealCO2mass = Calculator(registrationName="seal CO2 mass", Input=rktonnodes)
sealCO2mass.ResultArrayName = "SealTot [kg]"
sealCO2mass.Function = '("node CO2 mass fraction in gas"*"node gas mass density"*"node gas saturation"+"node CO2 mass fraction in liquid"*"node liquid mass density"*(1-"node gas saturation"))*(rocktype==1)'

# create a new 'Integrate Variables'
sealTot = IntegrateVariables(registrationName="SealTot", Input=sealCO2mass)

# create a new 'Plot Data Over Time'
dataovertime_2 = PlotDataOverTime(registrationName="data over time", Input=sealTot)
dataovertime_2.OnlyReportSelectionStatistics = 0

# create a new 'Calculator'
massCO2inboundary = Calculator(
    registrationName="mass CO2 in boundary (rocktypes 2-5)", Input=rktonnodes
)
massCO2inboundary.ResultArrayName = "BoundaryCO2 [kg]"
massCO2inboundary.Function = '("node CO2 mass fraction in gas"*"node gas mass density"*"node gas saturation"+"node CO2 mass fraction in liquid"*"node liquid mass density"*(1-"node gas saturation"))* ((rocktype >= 2) * (rocktype <= 5))'

# create a new 'Integrate Variables'
BoundaryCO2 = IntegrateVariables(
    registrationName="BoundaryCO2", Input=massCO2inboundary
)

# create a new 'Plot Data Over Time'
dataovertime_3 = PlotDataOverTime(registrationName="data over time2", Input=BoundaryCO2)
dataovertime_3.OnlyReportSelectionStatistics = 0

# create a new 'Gradient'
gradient1 = Gradient(registrationName="Gradient1", Input=stateswithmassandrocktypes)
gradient1.ScalarArray = ["POINTS", "node CO2 solubility ratio"]

# create a new 'Calculator'
calculator1 = Calculator(registrationName="Calculator1", Input=gradient1)
calculator1.ResultArrayName = "M_C [kg]"
calculator1.Function = "(Gradient[0]^2 + Gradient[1]^2 + Gradient[2]^2)^.5"

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(
    registrationName="IntegrateVariables1", Input=calculator1
)

# create a new 'Plot Data Over Time'
plotDataOverTime1 = PlotDataOverTime(
    registrationName="PlotDataOverTime1", Input=integrateVariables1
)
plotDataOverTime1.OnlyReportSelectionStatistics = 0

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
cSV1 = CreateExtractor("CSV", dataovertime_1, registrationName="CSV1")
# init the 'CSV' selected for 'Writer'
cSV1.Writer.FileName = "Box-A.csv"
cSV1.Writer.FieldAssociation = "Row Data"

# create extractor
cSV2 = CreateExtractor("CSV", dataovertime, registrationName="CSV2")
# init the 'CSV' selected for 'Writer'
cSV2.Writer.FileName = "Box-B.csv"
cSV2.Writer.FieldAssociation = "Row Data"

# create extractor
cSV4 = CreateExtractor("CSV", pressureovertime, registrationName="CSV4")
# init the 'CSV' selected for 'Writer'
cSV4.Writer.FileName = "POP2.csv"
cSV4.Writer.FieldAssociation = "Row Data"

# create extractor
cSV3 = CreateExtractor("CSV", pressureovertime_1, registrationName="CSV3")
# init the 'CSV' selected for 'Writer'
cSV3.Writer.FileName = "POP1.csv"
cSV3.Writer.FieldAssociation = "Row Data"

# create extractor
cSV5 = CreateExtractor("CSV", observation, registrationName="CSV5")
# init the 'CSV' selected for 'Writer'
cSV5.Writer.FileName = "dense_data_{time:6e}.csv"
# trace defaults for the extractor.
cSV5.Trigger = "Time Step"

# create extractor
cSV6 = CreateExtractor("CSV", dataovertime_2, registrationName="CSV1")
# init the 'CSV' selected for 'Writer'
cSV6.Writer.FileName = "SealTot.csv"
# trace defaults for the extractor.
cSV6.Trigger = "Time Step"
cSV6.Writer.FieldAssociation = "Row Data"

# create extractor
cSV7 = CreateExtractor("CSV", plotDataOverTime1, registrationName="CSV2")
# init the 'CSV' selected for 'Writer'
cSV7.Writer.FileName = "convC.csv"
# trace defaults for the extractor.
cSV7.Trigger = "Time Step"
cSV7.Writer.FieldAssociation = "Row Data"

# create extractor
cSV8 = CreateExtractor("CSV", dataovertime_3, registrationName="CSV3")
# init the 'CSV' selected for 'Writer'
cSV8.Writer.FileName = "BoundaryCO2.csv"
# trace defaults for the extractor.
cSV8.Trigger = "Time Step"
cSV8.Writer.FieldAssociation = "Row Data"

##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
SaveExtracts(ExtractsOutputDirectory=str(extract_directory))
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://kitware.github.io/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
