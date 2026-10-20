"""Running the auxiliary EPSR programs (setup and analysis tools) from Python.

Run from the ``python/`` folder:

    uv run python examples/04_tools/auxiliary_tools.py    # a few seconds

EPSR comes with 35 auxiliary programs besides the refinement itself: setup programs
(``makemole``, ``mixato``, ``nwts``...) and analysis programs (``partials``, ``coord``,
``sdf``...). They run on a copy of the DMSO/water test case (``template/``, from ``EPSR-test/``) in three
ways: keyword-file programs with a Python dict, programs with command-line arguments,
and programs that ask their questions on standard input. Everything goes to
``output/`` next to this script.
"""

import os
from pathlib import Path

import epsr
from epsr import tools

EXAMPLE_FOLDER = Path(__file__).resolve().parent
TEMPLATE = (
    EXAMPLE_FOLDER / "template"
)  # the run files of EPSR-test/: 200 DMSO + 400 water
OUTPUT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLE_FOLDER)) / "output"

project = epsr.Project.from_template(TEMPLATE, OUTPUT / "epsr_run", overwrite=True)
folder = project.folder  # the tools work on the folder, as EPSR's own programs do
print(len(tools.commands()), "tools:", ", ".join(tools.commands()))

# 1. A keyword-file setup tool: neutron weights from an existing .NWTS.dat file.
tools.nwts(folder, "DMSODWaterD.NWTS.dat")
print("nwts wrote", (folder / "DMSODWaterD.NWTStot.wts").name)

# 2. Write an analysis input file from a dict, then run the analysis.
partials_keywords = tools.keywords("partials")
print("partials keywords:", partials_keywords.global_keywords)
tools.write_keyword_file(
    folder / "water.PARTIALS.dat",
    "partials",
    {
        "fnameato": "200DMSO400Water.ato",
        "nz": 300,
        "atom-c": "O1",
        "atom-z": "*",
        "atom-r": "*",
    },
)
tools.partials(folder, "water.PARTIALS.dat")
partials_table = epsr.read_table(folder / "water.PARTIALS.g01")
print(
    f"partials around water O: {len(partials_table.labels)} columns, "
    f"r up to {partials_table.x[-1]:.1f} Å"
)
print("  first labels:", partials_table.labels[:4])

# 3. Tools that take command-line arguments or interactive answers.
#    writexyz converts an .ato configuration to .xyz. With arguments on the command line:
#    whole box, no rotation, nothing excluded, numbered output files.
whole_box_arguments = [
    "200DMSO400Water",
    "noappend",
    "0",
    "100",
    "100",
    "-100",
    "100",
    "0",
    "0",
    "0",
    "0",
]
tools.run(folder, "writexyz", whole_box_arguments)
print(
    "writexyz wrote", sorted(path.name for path in folder.glob("200DMSO400Water0*.xyz"))
)

#    The same program asking its questions on standard input: append? centre molecule,
#    x/y and z limits, Euler angles, atom types to exclude. One answer per line.
writexyz_answers = ["Y", "1", "10 10", "-10 10", "10 20 30", "3 4"]
tools.run(folder, "writexyz", ["DMSO.ato"], answers=writexyz_answers)
print("writexyz (interactive) wrote", (folder / "DMSOappend.xyz").name)

# 4. Errors are Python exceptions.
try:
    tools.nwts(folder, "does_not_exist.NWTS.dat")
except FileNotFoundError as error:
    print("error:", error)
