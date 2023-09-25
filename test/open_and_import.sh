#!/usr/bin/env bash

# given a path to a VSProject, open it, create a new AMP evaluation, import all variants that pass filters, then save the project

project_dir=$1
project_base=$(dirname "$project_dir")

vspipeline \
  -c project_open "$project_dir" \
  -c set_current_workflow workflow="amp" id="Workflow2" \
  -c run_workflow_script "/opt/VarSeq-2.3.0/Data/ReportTemplates/Test Import All Variants/script_cancer.js" "$project_base" \
  -c project_save \
  -c project_close
