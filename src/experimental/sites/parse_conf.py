import yaml

with open("conf.yaml") as f:
    conf = yaml.safe_load(f)
print(conf)

param = lambda name, value: f"integer, parameter :: {name} = {value}"
sites = [param(f"{name}_site", k + 1) for k, name in enumerate(conf["sites"])]
sites.append(param("nb_sites", len(conf["sites"])))
sites_declaration = "\n".join(sites)

with open("generated-sites.F90", "w") as f:
    f.write(
        f"""module Sites
implicit none
{sites_declaration}
end module Sites
"""
    )

param = lambda name, value: f"{name} = {value}"
sites = [param(f"{name}_site", k) for k, name in enumerate(conf["sites"])]
sites.append(param("nb_sites", len(conf["sites"])))
sites_declaration = ",\n".join(sites)

with open("generated-sites.h", "w") as f:
    f.write(
        f"""#pragma once
namespace ComPASS
{{
enum Sites: int {{
{sites_declaration}
}};
}} // namespace ComPASS
"""
    )
