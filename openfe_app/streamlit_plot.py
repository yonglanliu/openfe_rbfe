import matplotlib.pyplot as plt
from io import BytesIO
import numpy as np
import matplotlib as mpl

def plot_fr_convergence(df,
                        title="Forward and Reverse free energy convergence",
                        xlabel="Fraction of uncorrelated samples",
                        ylabel="ΔG (kcal/mol)",
                        show_errorbars=True,
                        marker_size=7,
                        line_width=2.5,
                        cap_size=3,
                        alpha_band=0.15,
                        show_band=True,
                        band_center=None,
                        band_halfwidth=0.05,
                        xlim=(0.0, 1.02),
                        ylim=None):
    fig, ax = plt.subplots(figsize=(7.5, 5))

    x = df["fraction"].values
    f = df["forward_DG (kcal/mol)"].values
    r = df["reverse_DG (kcal/mol)"].values
    fe = df["forward_error (kcal/mol)"].values
    re = df["reverse_error (kcal/mol)"].values

    if show_errorbars:
        ax.errorbar(x, f, yerr=fe, marker="o", color="red", markersize=marker_size, alpha=0.3,
                    linewidth=line_width, capsize=cap_size, label="Forward")
        ax.errorbar(x, r, yerr=re, marker="o", color="blue", markersize=marker_size, alpha=0.3,
                    linewidth=line_width, capsize=cap_size, label="Reverse")
    else:
        ax.plot(x, f, marker="o", alpha=0.9, markersize=marker_size, linewidth=line_width, label="Forward")
        ax.plot(x, r, marker="o", alpha=0.9, markersize=marker_size, linewidth=line_width, label="Reverse")

    # Optional horizontal band (like your purple band)
    if show_band:
        if band_center is None:
            # default: use last forward point as "final" reference
            band_center = float(f[-1])
        ax.axhspan(band_center - band_halfwidth,
                   band_center + band_halfwidth,
                   alpha=alpha_band,
                   color="blue")

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.legend()
    ax.grid(False)
    fig.tight_layout()

    return fig

def plot_replica_exchange_matrix(M, 
                                 title="Replica exchange transition matrix",
                                 annotate=True, 
                                 fmt="{:.2f}",
                                 vmin=0.0, vmax=None):
    """
    M: (K,K) transition matrix
    """
    K = M.shape[0]
    fig, ax = plt.subplots(figsize=(6, 6))

    # heatmap
    im = ax.imshow(M, origin="upper", cmap="Greys", vmin=vmin, vmax=vmax)

    # ticks / labels
    ax.set_xticks(range(K))
    ax.set_yticks(range(K))
    ax.set_xticklabels(range(K))
    ax.set_yticklabels(range(K))
    ax.set_title(title)
    ax.set_xlabel("λ")
    ax.set_ylabel("λ")

    # gridlines to match the style
    ax.set_xticks(np.arange(-.5, K, 1), minor=True)
    ax.set_yticks(np.arange(-.5, K, 1), minor=True)
    ax.grid(which="minor", linestyle="-", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    # annotate values
    if annotate:
        norm = mpl.colors.Normalize(
            vmin=vmin,
            vmax=np.nanmax(M) if vmax is None else vmax
        )

        for i in range(K):
            for j in range(K):
                val = M[i, j]
                if val == 0:
                    continue

                # Choose text color based on background intensity
                text_color = "white" if norm(val) > 0.5 else "black"

                ax.text(
                    j, i, fmt.format(val),
                    ha="center", va="center",
                    fontsize=9,
                    color=text_color
                )
    return fig

def plot_mbar_matrix(M, 
                    title="MBAR overlap matrix",
                    annotate=True, 
                    fmt="{:.2f}",
                    vmin=0.0, vmax=None):
    """
    M: (K,K) transition matrix
    """
    K = M.shape[0]
    fig, ax = plt.subplots(figsize=(6, 6))

    # heatmap
    im = ax.imshow(M, origin="upper", cmap="Greys", vmin=vmin, vmax=vmax)

    # ticks / labels
    ax.set_xticks(range(K))
    ax.set_yticks(range(K))
    ax.set_xticklabels(range(K))
    ax.set_yticklabels(range(K))
    ax.set_title(title)
    ax.set_xlabel("λ")
    ax.set_ylabel("λ")

    # gridlines to match the style
    ax.set_xticks(np.arange(-.5, K, 1), minor=True)
    ax.set_yticks(np.arange(-.5, K, 1), minor=True)
    ax.grid(which="minor", linestyle="-", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    # annotate values
    if annotate:
        norm = mpl.colors.Normalize(
            vmin=vmin,
            vmax=np.nanmax(M) if vmax is None else vmax
        )

        for i in range(K):
            for j in range(K):
                val = M[i, j]
                if val == 0:
                    continue

                # Choose text color based on background intensity
                text_color = "white" if norm(val) > 0.5 else "black"

                ax.text(
                    j, i, fmt.format(val),
                    ha="center", va="center",
                    fontsize=9,
                    color=text_color
                )
    return fig


def fig_to_png_bytes(fig, dpi=300):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()