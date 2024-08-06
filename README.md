\documentclass{article}
\usepackage{hyperref}

\begin{document}

\title{PRIMEloci}
\maketitle

\section*{Overview}

PRIMEloci is a collection of scripts designed for processing and analyzing CAGE-seq data, predicting transcription start site (TSS) profiles, and filtering prediction results. This repository provides a flexible and modular approach to run various steps in the analysis pipeline, allowing users to execute specific steps or the entire pipeline based on their needs.

\tableofcontents

\section{Installation}

Clone the repository:

\begin{verbatim}
git clone https://github.com/yourusername/PRIMEloci.git
cd PRIMEloci
\end{verbatim}

Make the main script executable:

\begin{verbatim}
chmod +x run_scripts.sh
\end{verbatim}

\section{Usage}

\subsection{Configuration}

All the parameter settings are stored in the \texttt{config.sh} file. Modify this file to match your data paths and desired parameters.

\subsection{Running Scripts}

To run the scripts, use the \texttt{run_scripts.sh} file. This script accepts options to specify which steps to run.

\subsubsection{Examples}

Run all steps:

\begin{verbatim}
./run_scripts.sh --all
\end{verbatim}

Run specific steps:

\begin{verbatim}
./run_scripts.sh -1 -3 -4
\end{verbatim}

\section{Scripts}

\subsection{get\_ctss\_from\_bw.R}

Extracts CAGE-seq data from bigWig files.

\textbf{Usage:}

\begin{verbatim}
Rscript _get_ctss_from_bw.r -i <CAGE_DIR> -m <DESIGN_MATRIX> -o <OUTPUT_DIR> -c <CTSS_RSE_NAME> -k
\end{verbatim}

\subsection{get\_tc\_grl.R}

Generates TSS cluster data from extracted CAGE-seq data.

\textbf{Usage:}

\begin{verbatim}
Rscript _get_tc_from_ctss.r -c <OUTPUT_DIR>/<CTSS_RSE_NAME> -o <OUTPUT_DIR> -t <TC_GRL_NAME> -e <EXTENSION_DISTANCE>
\end{verbatim}

\subsection{get\_tc\_profiles.R}

Profiles the TSS clusters.

\textbf{Usage:}

\begin{verbatim}
Rscript _get_tc_profiles.r -c <OUTPUT_DIR>/<CTSS_RSE_NAME> -t <OUTPUT_DIR>/<TC_GRL_NAME> -o <OUTPUT_DIR> -n <PROFILE_MAIN_DIR> -r <PROFILE_SUB_DIR>
\end{verbatim}

\subsection{predict\_profile\_probabilities.py}

Predicts TSS profile probabilities using a pre-trained model.

\textbf{Usage:}

\begin{verbatim}
python3 _predict_profile_probabilities.py -w <SCRIPT_DIR> -m <MODEL_PATH> -p <OUTPUT_DIR>/<PROFILE_MAIN_DIR> -r <PROFILE_SUB_DIR> -n <PREFIX_OUT_NAME> -t <THRESHOLD>
\end{verbatim}

\subsection{filter\_bed\_to\_reduced\_gr.R}

Filters prediction results to produce a reduced set of genomic ranges.

\textbf{Usage:}

\begin{verbatim}
Rscript _filter_bed_to_reduced_gr.r -i <FILE>
\end{verbatim}

\section{Contributing}

Contributions are welcome! Please submit a pull request or open an issue to discuss any changes or improvements.

\section{License}

This project is licensed under the MIT License - see the \href{LICENSE}{LICENSE} file for details.

\end{document}
