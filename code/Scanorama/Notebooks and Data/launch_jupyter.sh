#!/bin/sh
echo "->We are now running the analysis."
# echo "Now running jupyter nbconvert --to notebook --execute mynotebook.ipynb"
# jupyter nbconvert --to notebook --execute /Users/edjuaro/GoogleDrive/molecular_tumor_board/rady_on_docker/rady/rady_patient_rnaseq_pipeline.ipynb
# pwd
# ls -l /build/kallisto/kallisto_linux-v0.44.0
cd /Notebooks
# jupyter notebook --ip 0.0.0.0 --no-browser --allow-root ## Running a notebook
# jupyter nbextension enable collapsible_headings/main  ## Enabling an extension
# jupyter nbextension enable --user collapsible_headings/main

echo "You may have dropped this from your clipboard:"
echo "\t$pbpaste"
echo "USE THIS URL:"
echo "\t\thttp://127.0.0.1:8888/"
echo "http://127.0.0.1:8888/" | pbcopy
jupyter-notebook --no-browser --port 8888 --ip=0.0.0.0 --NotebookApp.password='' --NotebookApp.token='' --NotebookApp.password_required=False --allow-root
# jupyter notebook --ip 0.0.0.0 --allow-root ## Running a notebook
