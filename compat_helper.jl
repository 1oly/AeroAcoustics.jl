# Stand-alone script to manually check compat when not using GitHub workflow
using CompatHelper
deps = CompatHelper.get_project_deps("Project.toml")
CompatHelper.get_latest_version_from_registries!(deps,CompatHelper.DEFAULT_REGISTRIES)