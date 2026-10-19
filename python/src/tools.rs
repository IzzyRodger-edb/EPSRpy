//! The 35 auxiliary EPSR programs (`epsr_core::tools`) and their keyword tables.

use crate::convert::to_py;
use epsr_core::fortran::TextFile;
use epsr_core::rng::{Ran1, IM};
use epsr_core::symtab::{
    EPSR_SYMBOLS, MARK_IEREQFIRST, MARK_IEREQLAST, MARK_IINITREF, MARK_INCOEFFP,
    MARK_IQMINREF, MARK_IRESETREF, MARK_IRMINFAC, MARK_IROVERLAP, MARK_NSYMBOL,
    MARK_NSYMBOLDATA, MARK_NSYMBOLUPTODATA,
};
use epsr_core::symtab_tools;
use epsr_core::tools::{nwts, run_tool as core_run_tool, ToolContext, COMMANDS};
use epsr_core::{EpsrError, State};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::path::PathBuf;

/// Run one auxiliary program in `folder`, as `epsr-rs <folder> <command> <args...>`.
///
/// `answers` is the text the Fortran program would read from standard input.
/// Files in `folder` are created or rewritten in place. Releases the GIL.
#[pyfunction]
#[pyo3(signature = (folder, command, args=Vec::new(), answers=String::new(), verbose=false))]
pub fn run_tool(
    py: Python<'_>,
    folder: PathBuf,
    command: &str,
    args: Vec<String>,
    answers: String,
    verbose: bool,
) -> PyResult<()> {
    let command = command.to_ascii_lowercase();
    if !COMMANDS.contains(&command.as_str()) {
        return Err(PyValueError::new_err(format!(
            "unknown EPSR tool {command:?}; valid tools: {}",
            COMMANDS.join(", ")
        )));
    }
    py.detach(move || {
        let standard_input = TextFile::from_string(&answers);
        let mut context = ToolContext::new(&folder, args, standard_input, verbose);
        core_run_tool(&command, &mut context)
    })
    .map_err(to_py)
}

/// Names of the auxiliary programs accepted by `run_tool`.
#[pyfunction]
pub fn tool_commands() -> Vec<&'static str> {
    COMMANDS.to_vec()
}

/// Names of the programs that have a keyword-file table (`tool_keywords`).
#[pyfunction]
pub fn tool_tables() -> Vec<&'static str> {
    symtab_tools::ALL.iter().map(|table| table.name).collect()
}

/// The keyword table of a keyword-file driven program (e.g. `nwts`, `partials`).
///
/// Returns a dict with `name`, `extinp`, `prompt`, `dataprompt`,
/// `dataprompttext`, `symbols` (list of `(symbol, type, iflag, default,
/// description)`) and `marks` (dict, e.g. `nsymboluptodata`: number of
/// keywords before the first repeated `dataprompt` block).
#[pyfunction]
pub fn tool_keywords<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyDict>> {
    let table = symtab_tools::ALL
        .iter()
        .find(|table| table.name.eq_ignore_ascii_case(name))
        .ok_or_else(|| {
            PyValueError::new_err(format!(
                "no keyword table for {name:?}; tables exist for: {}",
                tool_tables().join(", ")
            ))
        })?;
    let description = PyDict::new(py);
    description.set_item("name", table.name)?;
    description.set_item("extinp", table.extinp)?;
    description.set_item("prompt", table.prompt)?;
    description.set_item("dataprompt", table.dataprompt)?;
    description.set_item("dataprompttext", table.dataprompttext)?;
    let symbols = PyList::empty(py);
    for (symbol, value_type, flag, default, symbol_description) in table.symbols {
        symbols.append((
            *symbol,
            value_type.to_string(),
            *flag,
            *default,
            *symbol_description,
        ))?;
    }
    description.set_item("symbols", symbols)?;
    let marks = PyDict::new(py);
    for (mark_name, position) in table.marks {
        marks.set_item(*mark_name, *position)?;
    }
    description.set_item("marks", marks)?;
    Ok(description)
}

/// The file an `.EPSR.inp` keyword is saved to, from its 1-based table position.
///
/// Mirrors the engine's writers (`save_inp_file` and `save_pcof_file` in
/// `symbols.rs`): the `.pcof` header holds positions `iereqfirst+1..=iereqlast`
/// and `iqminref+1..nsymboluptodata`, the `data n` blocks the positions after
/// `nsymboluptodata`, and the `.inp` header everything else.
fn epsr_keyword_file(position: i32) -> &'static str {
    if position > MARK_NSYMBOLUPTODATA {
        "data"
    } else if (position > MARK_IEREQFIRST && position <= MARK_IEREQLAST)
        || (position > MARK_IQMINREF && position < MARK_NSYMBOLUPTODATA)
    {
        "pcof"
    } else {
        "inp"
    }
}

/// The `.EPSR.inp` keyword table (`symtab.rs`, `EPSR_SYMBOLS`), which also holds
/// the `.pcof` header keywords.
///
/// Returns a dict with `symbols` (list of `(symbol, type, iflag, default,
/// description, file)`, where `file` is `"inp"`, `"pcof"` or `"data"`, the
/// file or block the engine saves the keyword to) and `marks` (the table
/// positions of `symtab.rs`, e.g. `nsymboluptodata`).
#[pyfunction]
pub fn epsr_keywords<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let description = PyDict::new(py);
    let symbols = PyList::empty(py);
    for (index, (symbol, value_type, flag, default, symbol_description)) in
        EPSR_SYMBOLS.iter().enumerate()
    {
        symbols.append((
            *symbol,
            value_type.to_string(),
            *flag,
            *default,
            *symbol_description,
            epsr_keyword_file(index as i32 + 1),
        ))?;
    }
    description.set_item("symbols", symbols)?;
    let marks = PyDict::new(py);
    for (mark_name, position) in [
        ("iereqfirst", MARK_IEREQFIRST),
        ("iereqlast", MARK_IEREQLAST),
        ("iresetref", MARK_IRESETREF),
        ("iinitref", MARK_IINITREF),
        ("iqminref", MARK_IQMINREF),
        ("iroverlap", MARK_IROVERLAP),
        ("irminfac", MARK_IRMINFAC),
        ("incoeffp", MARK_INCOEFFP),
        ("nsymboluptodata", MARK_NSYMBOLUPTODATA),
        ("nsymboldata", MARK_NSYMBOLDATA),
        ("nsymbol", MARK_NSYMBOL),
    ] {
        marks.set_item(mark_name, position)?;
    }
    description.set_item("marks", marks)?;
    Ok(description)
}

/// The neutron scattering-length table used by `nwts` (Sears 1991):
/// a list of `(symbol, mass_number, scattering_length)`, mass number 0 being
/// the natural element, lengths in 10⁻¹² cm.
#[pyfunction]
pub fn isotope_table() -> Vec<(&'static str, i32, f32)> {
    nwts::isotope_table().to_vec()
}

/// Seed the random-number generator of a box `.ato` file.
///
/// Loads the file with the engine's reader (`minitc`), sets the generator
/// state to `idum = -seed` with every other entry zero (the standard way to
/// initialise `ran1`), and writes the file back with the engine's writer
/// (`moloutc`), as `randomise` does. Valid seeds are 1 to 2³¹ − 2; the
/// all-zero state that `mixato` leaves behaves like seed 1. Releases the GIL.
#[pyfunction]
pub fn set_rng_seed(py: Python<'_>, ato_path: PathBuf, seed: i64) -> PyResult<()> {
    let largest_seed = i64::from(IM) - 1;
    if !(1..=largest_seed).contains(&seed) {
        return Err(PyValueError::new_err(format!(
            "seed must be between 1 and {largest_seed}, got {seed}"
        )));
    }
    py.detach(move || {
        let file_name = ato_path.to_string_lossy().into_owned();
        let mut state = State::new();
        state.mol.ntype = 0;
        state.mol.nmoltype = 0;
        state.mol.sizefactor = 1.0;
        state.minitc(&file_name)?;
        if state.dist.rho == 0.0 {
            return Err(EpsrError::Input(format!(
                "{file_name} is not a box .ato file (zero density)"
            )));
        }
        state.rng = Ran1 {
            idum: -(seed as i32),
            ..Ran1::default()
        };
        state.moloutc(&file_name)
    })
    .map_err(to_py)
}
