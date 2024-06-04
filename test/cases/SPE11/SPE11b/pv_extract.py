""" python pv_extract.py -h
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

# create a new 'Plane'
boxCgeometry = Plane(registrationName="Box C geometry")
boxCgeometry.Origin = [3300.0, 50.0, 100.0]
boxCgeometry.Point1 = [7800.0, 50.0, 100.0]
boxCgeometry.Point2 = [3300.0, 50.0, 400.0]
boxCgeometry.XResolution = 900
boxCgeometry.YResolution = 60

# create a new 'Plane'
boxBgeometry = Plane(registrationName="Box B geometry")
boxBgeometry.Origin = [100.0, 50.0, 600.0]
boxBgeometry.Point1 = [3300.0, 50.0, 600.0]
boxBgeometry.Point2 = [100.0, 50.0, 1200.0]
boxBgeometry.XResolution = 640
boxBgeometry.YResolution = 120

# create a new 'Plane'
observationgeometry = Plane(registrationName="observation geometry")
observationgeometry.Origin = [0.0, 50.0, 0.0]
observationgeometry.Point1 = [8400.0, 50.0, 0.0]
observationgeometry.Point2 = [0.0, 50.0, 1200.0]
observationgeometry.XResolution = 840
observationgeometry.YResolution = 120

# create a new 'XML Unstructured Grid Reader'
rocktypes = XMLUnstructuredGridReader(
    registrationName="rocktypes", FileName=[str(rocktype_vtu_file)]
)
rocktypes.CellArrayStatus = ["rkt"]
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
massCO2ingas.Function = '"node CO2 mass fraction in gas"*"node gas mass density"'

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
massCO2inliquid.Function = (
    '"node CO2 mass fraction in liquid"*"node liquid mass density"'
)

# create a new 'Calculator'
massCO2insealrtk1 = Calculator(
    registrationName="mass CO2 in seal (rtk==1)", Input=massCO2inliquid
)
massCO2insealrtk1.ResultArrayName = "sealB"
massCO2insealrtk1.Function = '("node CO2 mass fraction in liquid"*"node liquid mass density"+"node CO2 mass fraction in gas"*"node gas mass density")*(rkt==1)'

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
pOP2.ProbeType.Center = [5100.0, 50.0, 1100.0]

# create a new 'Plot Data Over Time'
pressureovertime = PlotDataOverTime(registrationName="pressure over time", Input=pOP2)
pressureovertime.OnlyReportSelectionStatistics = 0

# create a new 'Plane'
boxAgeometry = Plane(registrationName="Box A geometry")
boxAgeometry.Origin = [3300.0, 50.0, 0.0]
boxAgeometry.Point1 = [8300.0, 50.0, 0.0]
boxAgeometry.Point2 = [3300.0, 50.0, 600.0]
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
pOP1.ProbeType.Center = [4500.0, 50.0, 500.0]

# create a new 'Plot Data Over Time'
pressureovertime_1 = PlotDataOverTime(registrationName="pressure over time", Input=pOP1)
pressureovertime_1.OnlyReportSelectionStatistics = 0

# create a new 'Calculator'
massC02ingas = Calculator(registrationName="mass C02 in gas", Input=boxA)
massC02ingas.ResultArrayName = "mCO2 gas"
massC02ingas.Function = '"node CO2 mass fraction in gas"*"node gas mass density"'

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
massCO2inliquid_1.Function = (
    '"node CO2 mass fraction in liquid"*"node liquid mass density"'
)

# create a new 'Calculator'
massCO2insealrtk1_1 = Calculator(
    registrationName="mass CO2 in seal (rtk==1)", Input=massCO2inliquid_1
)
massCO2insealrtk1_1.ResultArrayName = "sealA"
massCO2insealrtk1_1.Function = '("node CO2 mass fraction in liquid"*"node liquid mass density"+"node CO2 mass fraction in gas"*"node gas mass density")*(rkt==1)'

# create a new 'Integrate Variables'
integration_1 = IntegrateVariables(
    registrationName="integration", Input=massCO2insealrtk1_1
)

# create a new 'Plot Data Over Time'
dataovertime_1 = PlotDataOverTime(
    registrationName="data over time", Input=integration_1
)
dataovertime_1.OnlyReportSelectionStatistics = 0

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
cSV1 = CreateExtractor("CSV", dataovertime_1, registrationName="CSV1")
# trace defaults for the extractor.
cSV1.Trigger = "Time Step"

# init the 'CSV' selected for 'Writer'
cSV1.Writer.FileName = "Box-A.csv"
cSV1.Writer.FieldAssociation = "Row Data"

# create extractor
cSV2 = CreateExtractor("CSV", dataovertime, registrationName="CSV2")
# trace defaults for the extractor.
cSV2.Trigger = "Time Step"

# init the 'CSV' selected for 'Writer'
cSV2.Writer.FileName = "Box-B.csv"
cSV2.Writer.FieldAssociation = "Row Data"

# create extractor
cSV4 = CreateExtractor("CSV", pressureovertime, registrationName="CSV4")
# trace defaults for the extractor.
cSV4.Trigger = "Time Step"

# init the 'CSV' selected for 'Writer'
cSV4.Writer.FileName = "POP2.csv"
cSV4.Writer.FieldAssociation = "Row Data"

# create extractor
cSV3 = CreateExtractor("CSV", pressureovertime_1, registrationName="CSV3")
# trace defaults for the extractor.
cSV3.Trigger = "Time Step"

# init the 'CSV' selected for 'Writer'
cSV3.Writer.FileName = "POP1.csv"
cSV3.Writer.FieldAssociation = "Row Data"

# ----------------------------------------------------------------
# restore active source
SetActiveSource(observation)
# ----------------------------------------------------------------


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
