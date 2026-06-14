from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st


def read_dg_results(tsv_file, lig_col, dg_col, unc_col):
    df = pd.read_csv(str(tsv_file), sep="\t")

    required_cols = [lig_col, dg_col, unc_col]
    missing_cols = [c for c in required_cols if c not in df.columns]

    if missing_cols:
        raise ValueError(
            f"Missing columns: {missing_cols}\n"
            f"Available columns: {list(df.columns)}"
        )

    df[dg_col] = pd.to_numeric(df[dg_col], errors="coerce")
    df[unc_col] = pd.to_numeric(df[unc_col], errors="coerce")

    error_rows = df[df[[dg_col, unc_col]].isna().any(axis=1)].copy()

    plot_df = df.dropna(subset=[dg_col, unc_col]).copy()
    plot_df = plot_df.sort_values(dg_col).reset_index(drop=True)

    return df, plot_df, error_rows


def plot_dg_mle(
    plot_df,
    lig_col,
    dg_col,
    unc_col,
    fig_width=12,
    fig_height=6,
    marker_size=7,
    capsize=5,
    line_width=1.5,
    label_rotation=45,
    x_font_size=14,
    y_font_size=14,
    ylabel_font_size=20,
    point_color="black",
    zero_line_color="red",
    title="",
):
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    ax.errorbar(
        plot_df[lig_col],
        plot_df[dg_col],
        yerr=plot_df[unc_col],
        fmt="o",
        markersize=marker_size,
        capsize=capsize,
        linewidth=line_width,
        color=point_color,
        alpha=0.8,
    )

    ax.axhline(0, linestyle="--", color=zero_line_color)

    ax.set_ylabel(r"MLE $\Delta G$ (kcal/mol)", fontsize=ylabel_font_size)

    if title:
        ax.set_title(title)

    ax.tick_params(axis="x", labelsize=x_font_size)
    ax.tick_params(axis="y", labelsize=y_font_size)

    plt.setp(
        ax.get_xticklabels(),
        rotation=label_rotation,
        ha="right",
    )

    fig.tight_layout()

    return fig


def tab3_design(analysis_tsv_picker, workdir):
    st.header("DG Trend")

    selected_tsv_file = analysis_tsv_picker(
        label="Input DG results TSV file",
        start_dir=str(workdir),
        allowed_extensions=(".tsv",),
        key_prefix="dg_trend_tsv",
    )

    st.divider()

    if selected_tsv_file is None:
        st.info("Click **List**, browse to your dg_results.tsv file, then click **Use TSV file**.")
        return

    tsv_file = Path(selected_tsv_file)

    if not tsv_file.exists():
        st.error(f"Cannot find TSV file:\n\n{tsv_file}")
        return

    left_col, right_col = st.columns([1, 3], gap="large")

    with left_col:
        st.subheader("Plot controls")

        # with st.expander("Column names", expanded=True):
        #     lig_col = st.text_input(
        #         "Ligand column",
        #         value="ligand",
        #         key="dg_lig_col",
        #     )

        #     dg_col = st.text_input(
        #         "DG column",
        #         value="DG(MLE) (kcal/mol)",
        #         key="dg_mle_col",
        #     )

        #     unc_col = st.text_input(
        #         "Uncertainty column",
        #         value="uncertainty (kcal/mol)",
        #         key="dg_unc_col",
        #     )


        lig_col = "ligand"
        dg_col ="DG(MLE) (kcal/mol)"
        unc_col = "uncertainty (kcal/mol)"

        with st.expander("Figure style", expanded=True):
            cc1, cc2 = st.columns(2)
            with cc1:
                fig_width = st.slider(
                    "Figure width",
                    min_value=6,
                    max_value=24,
                    value=12,
                    step=1,
                    key="dg_fig_width",
                )
            with cc2:
                fig_height = st.slider(
                    "Figure height",
                    min_value=4,
                    max_value=16,
                    value=6,
                    step=1,
                    key="dg_fig_height",
                )
            cd1, cd2 = st.columns(2)
            with cd1:
                marker_size = st.slider(
                    "Marker size",
                    min_value=3,
                    max_value=20,
                    value=7,
                    step=1,
                    key="dg_marker_size",
                )
            with cd2:
                capsize = st.slider(
                    "Error bar cap size",
                    min_value=0,
                    max_value=15,
                    value=5,
                    step=1,
                    key="dg_capsize",
                )
            ce1, ce2 = st.columns(2)
            with ce1:
                line_width = st.slider(
                    "Error bar line width",
                    min_value=0.5,
                    max_value=5.0,
                    value=1.5,
                    step=0.1,
                    key="dg_line_width",
                )
            with ce2:
                label_rotation = st.slider(
                    "X label rotation",
                    min_value=0,
                    max_value=90,
                    value=45,
                    step=5,
                    key="dg_label_rotation",
                )

            x_font_size = st.slider(
                "X label font size",
                min_value=6,
                max_value=30,
                value=14,
                step=1,
                key="dg_x_font_size",
            )

            cf1, cf2 = st.columns(2)
            with cf1:
                y_font_size = st.slider(
                    "Y tick font size",
                    min_value=6,
                    max_value=30,
                    value=14,
                    step=1,
                    key="dg_y_font_size",
                )
            with cf2:
                ylabel_font_size = st.slider(
                    "Y axis title font size",
                    min_value=8,
                    max_value=40,
                    value=20,
                    step=1,
                    key="dg_ylabel_font_size",
                )
            ch1, ch2 = st.columns(2)
            with ch1:
                point_color = st.color_picker(
                    "Point / error bar color",
                    value="#000000",
                    key="dg_point_color",
                )
            with ch2:
                zero_line_color = st.color_picker(
                    "Zero line color",
                    value="#FF0000",
                    key="dg_zero_line_color",
                )

        with st.expander("Output", expanded=True):
            if tsv_file.parent.name == "analysis":
                analysis_dir = tsv_file.parent
            else:
                analysis_dir = tsv_file.parent / "analysis"

            out_png = Path(
                st.text_input(
                    "Output PNG",
                    value=str(analysis_dir / "DG_MLE_plot.png"),
                    key="dg_out_png",
                )
            )

            out_pdf = Path(
                st.text_input(
                    "Output PDF",
                    value=str(analysis_dir / "DG_MLE_plot.pdf"),
                    key="dg_out_pdf",
                )
            )

            dpi = st.slider(
                "Save DPI",
                min_value=100,
                max_value=900,
                value=600,
                step=50,
                key="dg_dpi",
            )

            save_pdf = st.checkbox(
                "Also save PDF",
                value=False,
                key="dg_save_pdf",
            )

    try:
        df, plot_df, error_rows = read_dg_results(
            tsv_file=tsv_file,
            lig_col=lig_col,
            dg_col=dg_col,
            unc_col=unc_col,
        )
    except Exception as e:
        with right_col:
            st.error(f"Could not read DG results TSV:\n\n{e}")
        return

    if plot_df.empty:
        with right_col:
            st.error("No valid numeric DG rows found after removing Error/NaN rows.")
        return

    fig = plot_dg_mle(
        plot_df=plot_df,
        lig_col=lig_col,
        dg_col=dg_col,
        unc_col=unc_col,
        fig_width=fig_width,
        fig_height=fig_height,
        marker_size=marker_size,
        capsize=capsize,
        line_width=line_width,
        label_rotation=label_rotation,
        x_font_size=x_font_size,
        y_font_size=y_font_size,
        ylabel_font_size=ylabel_font_size,
        point_color=point_color,
        zero_line_color=zero_line_color,
    )

    with right_col:
        st.subheader("DG MLE plot")

        b1, b2, b3 = st.columns(3)

        with b1:
            if st.button("Save PNG", type="primary", key="dg_save_png_button"):
                try:
                    out_png.parent.mkdir(parents=True, exist_ok=True)
                    fig.savefig(str(out_png), dpi=dpi, bbox_inches="tight")
                    st.success(f"Saved PNG to:\n\n{out_png}")
                except Exception as e:
                    st.error(f"Could not save PNG:\n\n{e}")

        with b2:
            if st.button("Save PDF", type="primary", key="dg_save_pdf_button"):
                try:
                    out_pdf.parent.mkdir(parents=True, exist_ok=True)
                    fig.savefig(str(out_pdf), bbox_inches="tight")
                    st.success(f"Saved PDF to:\n\n{out_pdf}")
                except Exception as e:
                    st.error(f"Could not save PDF:\n\n{e}")

        with b3:
            st.info(f"Valid rows: {len(plot_df)} | Skipped rows: {len(error_rows)}")

        st.pyplot(fig)
        plt.close(fig)

    # with st.expander("Valid plotted rows", expanded=False):
    #     st.dataframe(plot_df, use_container_width=True)

    # if not error_rows.empty:
    #     with st.expander("Skipped Error / NaN rows", expanded=False):
    #         st.dataframe(error_rows, use_container_width=True)
