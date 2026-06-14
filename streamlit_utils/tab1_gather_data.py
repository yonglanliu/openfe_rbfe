import streamlit as st
import subprocess
from pathlib import Path

def run_command_live(cmd, success_message):
    """
    Run a subprocess command and show stdout/stderr live in Streamlit.
    """
    st.write("Running command:")
    st.code(" ".join(cmd))

    log_box = st.empty()
    logs = []

    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        if process.stdout is not None:
            for line in process.stdout:
                logs.append(line)
                log_box.code("".join(logs[-300:]))

        return_code = process.wait()

        if return_code == 0:
            st.success(success_message)
            return True

        st.error(f"Command failed with return code {return_code}")
        return False

    except FileNotFoundError:
        st.error(
            "Could not find `openfe`. Make sure Streamlit is running in the same "
            "conda environment where OpenFE is installed."
        )
        return False

    except Exception as e:
        st.error(f"Unexpected error:\n\n{e}")
        return False


def openfe_gather_panel(
    title,
    result_dir,
    analysis_dir,
    report_type,
    default_filename,
    text_key,
    button_key,
):
    """
    Reusable panel for running:

        openfe gather RESULT_DIR --report REPORT_TYPE -o OUTPUT_FILE
    """
    st.subheader(title)

    if result_dir is None:
        st.warning("Please choose the results directory first.")
        return

    if analysis_dir is None:
        st.warning("Please choose a directory to save data first.")
        return

    result_dir = Path(result_dir)
    analysis_dir = Path(analysis_dir)

    if not result_dir.exists():
        st.error(f"Results directory does not exist:\n\n{result_dir}")
        return

    analysis_dir.mkdir(parents=True, exist_ok=True)

    default_output_file = analysis_dir / default_filename

    last_analysis_key = f"{text_key}_last_analysis_dir"

    if (
        last_analysis_key not in st.session_state
        or st.session_state[last_analysis_key] != str(analysis_dir)
    ):
        st.session_state[text_key] = str(default_output_file)
        st.session_state[last_analysis_key] = str(analysis_dir)

    output_file = st.text_input(
        "Output file",
        key=text_key,
    )

    output_file = Path(output_file)

    st.write("Output will be saved to:")
    st.code(str(output_file))

    if st.button(title, type="primary", key=button_key):
        output_file.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "openfe",
            "gather",
            str(result_dir),
            "--report",
            report_type,
            "-o",
            str(output_file),
        ]
        if report_type == "dg":
            cmd.append("--allow-partial")

        run_command_live(
            cmd,
            success_message=f"Finished {report_type} data gathering:\n\n{output_file}",
        )



def tab1_design(directory_picker, workdir):

    st.header("Gather Data")

    st.info(
        "Before running data gathering, choose the OpenFE results directory "
        "and the directory where you want to save gathered TSV files."
    )

    result_dir = directory_picker(
        label="Custom Results Directory",
        start_dir=str(workdir),
        key_prefix="result_dir",
    )

    st.divider()

    analysis_dir = directory_picker(
        label="Directory to Save Data",
        start_dir=str(workdir),
        key_prefix="analysis_dir",
    )

    st.divider()


    c1, c2 = st.columns(2)

    with c1:
        openfe_gather_panel(
            title="Gather raw data",
            result_dir=result_dir,
            analysis_dir=analysis_dir,
            report_type="raw",
            default_filename="raw_results.tsv",
            text_key="raw_results_output_file",
            button_key="btn_gather_raw_data",
        )

    with c2:
        openfe_gather_panel(
            title="Gather dG data",
            result_dir=result_dir,
            analysis_dir=analysis_dir,
            report_type="dg",
            default_filename="dg_results.tsv",
            text_key="dg_results_output_file",
            button_key="btn_gather_dg_data",
        )
