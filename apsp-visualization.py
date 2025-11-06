import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path

st.set_page_config(
    page_title="APSP Performance Dashboard",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded",
)


@st.cache_data
def load_results(results_dir="apsp_results/results"):
    results_path = Path(results_dir)

    if not results_path.exists():
        return pd.DataFrame(), [f"Results directory not found: {results_dir}"]

    all_data = []
    errors = []

    for csv_file in results_path.glob("*.csv"):
        if csv_file.name.startswith("._"):
            continue

        try:
            df = pd.read_csv(csv_file)

            if df.empty or not {"algorithm", "matrix_size", "time_sec"}.issubset(
                df.columns
            ):
                errors.append(f"{csv_file.name}: Invalid file")
                continue

            # Handle processes column
            if "processes" not in df.columns:
                df["processes"] = df["algorithm"].apply(
                    lambda x: 1 if x in ["fw", "rs", "blocked"] else None
                )

            df[["time_sec", "matrix_size"]] = df[["time_sec", "matrix_size"]].apply(
                pd.to_numeric, errors="coerce"
            )
            if "processes" in df.columns:
                df["processes"] = pd.to_numeric(df["processes"], errors="coerce")

            df = df[
                (df["time_sec"] > 0)
                & df["time_sec"].notna()
                & df["matrix_size"].notna()
            ]

            if not df.empty:
                all_data.append(df)

        except Exception as e:
            errors.append(f"{csv_file.name}: {str(e)}")

    if not all_data:
        return pd.DataFrame(), errors + ["No valid data found"]

    return pd.concat(all_data, ignore_index=True), errors


def format_algorithm_name(alg):
    names = {
        "rs": "Repeated Squaring",
        "fw": "Floyd-Warshall (Sequential Baseline)",
        "blocked": "Blocked FW (Sequential)",
        "mpi": "MPI Blocking",
        "mpi-nb": "MPI Non-Blocking",
    }
    return names.get(alg, alg.upper())


def calculate_speedup(df):
    """Calculate speedup using sequential FW as baseline"""
    baseline = df[df["algorithm"] == "fw"]

    if baseline.empty:
        st.warning("No sequential FW baseline found - speedup calculations unavailable")
        return df.assign(speedup=1.0, efficiency=100.0)

    baseline_times = baseline.groupby("matrix_size")["time_sec"].median().reset_index()
    baseline_times.columns = ["matrix_size", "baseline_time"]

    result = df.merge(baseline_times, on="matrix_size", how="left")
    result["speedup"] = result["baseline_time"] / result["time_sec"]

    # Calculate efficiency only for parallel algorithms
    result["efficiency"] = result.apply(
        lambda row: (row["speedup"] / row["processes"]) * 100
        if pd.notna(row["processes"]) and row["processes"] > 1
        else None,
        axis=1,
    )

    return result


def main():
    st.title("APSP Performance Analysis")
    st.markdown("Performance analysis using **Sequential Floyd-Warshall** as baseline")

    df, errors = load_results()

    if df.empty:
        st.error("No benchmark data found")
        if errors:
            with st.expander("Show errors"):
                for error in errors:
                    st.text(error)
        st.info(
            "Place CSV files in the `results/` directory with columns: algorithm, matrix_size, time_sec, processes"
        )
        return

    st.success(f"Loaded {len(df)} benchmark runs")

    with st.expander("Data Summary"):
        st.write(f"**Algorithms found:** {', '.join(sorted(df['algorithm'].unique()))}")
        st.write(
            f"**Matrix sizes:** {', '.join(map(str, sorted(df['matrix_size'].unique())))}"
        )
        if "processes" in df.columns and df["processes"].notna().any():
            valid_procs = df[df["processes"].notna()]["processes"].unique()
            st.write(f"**Process counts:** {', '.join(map(str, sorted(valid_procs)))}")

        st.write("**Runs per algorithm:**")
        for algo, count in (
            df.groupby("algorithm").size().sort_values(ascending=False).items()
        ):
            st.write(f"  - {format_algorithm_name(algo)}: {count} runs")

    if errors:
        with st.expander(f"{len(errors)} warnings during loading"):
            for error in errors:
                st.text(error)

    st.sidebar.header("Filters")

    selected_algorithms = st.sidebar.multiselect(
        "Algorithms",
        options=sorted(df["algorithm"].unique()),
        default=sorted(df["algorithm"].unique()),
        format_func=format_algorithm_name,
    )

    selected_sizes = st.sidebar.multiselect(
        "Matrix Sizes",
        options=sorted(df["matrix_size"].unique()),
        default=sorted(df["matrix_size"].unique()),
    )

    filtered = df[
        df["algorithm"].isin(selected_algorithms)
        & df["matrix_size"].isin(selected_sizes)
    ]

    if filtered.empty:
        st.warning("No data matches filters")
        return

    st.header("Key Results")

    fastest = filtered.loc[filtered["time_sec"].idxmin()]
    speedup_df = calculate_speedup(filtered)

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Fastest Run", f"{fastest['time_sec']:.3f}s")
    col2.metric("Best Algorithm", format_algorithm_name(fastest["algorithm"]))

    proc_str = (
        f"p={int(fastest['processes'])}"
        if pd.notna(fastest.get("processes"))
        else "sequential"
    )
    col3.metric("Best Config", f"n={int(fastest['matrix_size'])}, {proc_str}")
    col4.metric("Max Speedup", f"{speedup_df['speedup'].max():.2f}×")

    st.markdown("---")

    st.header("1. Execution Time by Algorithm")

    summary = (
        filtered.groupby(["algorithm", "matrix_size"])
        .agg({"time_sec": "median"})
        .reset_index()
    )
    summary["algorithm_name"] = summary["algorithm"].apply(format_algorithm_name)

    fig1 = px.line(
        summary,
        x="matrix_size",
        y="time_sec",
        color="algorithm_name",
        markers=True,
        title="Median Execution Time (Sequential FW = Baseline)",
        labels={
            "time_sec": "Time (seconds)",
            "matrix_size": "Matrix Size (n)",
            "algorithm_name": "Algorithm",
        },
        height=450,
        log_y=True,
    )
    fig1.update_layout(
        legend=dict(
            orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, title=None
        ),
        hovermode="x unified",
    )
    st.plotly_chart(fig1, use_container_width=True)

    st.subheader("Performance Summary")
    perf_table = summary.pivot_table(
        values="time_sec", index="algorithm_name", columns="matrix_size"
    )

    def format_time(val):
        if pd.isna(val):
            return "-"
        return f"{val * 1000:.1f}ms" if val < 1 else f"{val:.3f}s"

    st.dataframe(
        perf_table.style.format(format_time).background_gradient(
            cmap="RdYlGn_r", axis=None
        ),
        use_container_width=True,
    )

    st.markdown("---")

    st.header("2. Speedup & Efficiency (vs Sequential FW)")

    # Filter for parallel algorithms only
    parallel_data = speedup_df[speedup_df["algorithm"].isin(["mpi", "mpi-nb"])].copy()

    if not parallel_data.empty and "processes" in parallel_data.columns:
        speedup_summary = (
            parallel_data.groupby(["algorithm", "matrix_size", "processes"])
            .agg({"speedup": "median", "efficiency": "median"})
            .reset_index()
        )
        speedup_summary["algorithm_name"] = speedup_summary["algorithm"].apply(
            format_algorithm_name
        )

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Speedup")
            speedup_table = speedup_summary.pivot_table(
                values="speedup",
                index=["algorithm_name", "matrix_size"],
                columns="processes",
            ).round(2)
            st.dataframe(
                speedup_table.style.background_gradient(cmap="Greens", axis=None),
                use_container_width=True,
            )

        with col2:
            st.subheader("Efficiency (%)")
            efficiency_table = speedup_summary.pivot_table(
                values="efficiency",
                index=["algorithm_name", "matrix_size"],
                columns="processes",
            ).round(1)
            st.dataframe(
                efficiency_table.style.background_gradient(cmap="Blues", axis=None),
                use_container_width=True,
            )
    else:
        st.info(
            "Speedup analysis requires sequential FW baseline and parallel algorithms with process counts"
        )

    st.markdown("---")

    st.header("3. Export Data")

    export_df = filtered.copy()
    export_df["algorithm"] = export_df["algorithm"].apply(format_algorithm_name)

    groupby_cols = ["algorithm", "matrix_size"]
    if "processes" in export_df.columns:
        groupby_cols.append("processes")

    export_summary = (
        export_df.groupby(groupby_cols)
        .agg({"time_sec": ["median", "mean", "std", "count"]})
        .reset_index()
    )

    col_names = ["Algorithm", "Matrix Size"]
    if "processes" in export_df.columns:
        col_names.append("Processes")
    col_names.extend(["Median Time", "Mean Time", "Std Dev", "Runs"])
    export_summary.columns = col_names

    st.download_button(
        label="Download Summary CSV",
        data=export_summary.to_csv(index=False),
        file_name="apsp_summary.csv",
        mime="text/csv",
    )

    st.dataframe(export_summary, use_container_width=True, height=400)


if __name__ == "__main__":
    main()
