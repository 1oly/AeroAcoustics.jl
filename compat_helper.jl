# Stand-alone script to manually check compat when not using GitHub workflow
using CompatHelper
dep_to_current_compat_entry,
    dep_to_current_compat_entry_verbatim,
    dep_to_latest_version,
    deps_with_missing_compat_entry = CompatHelper.get_project_deps("Project.toml")

CompatHelper.get_latest_version_from_registries!(dep_to_latest_version,CompatHelper.default_registries)