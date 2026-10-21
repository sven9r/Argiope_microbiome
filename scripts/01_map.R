# Figure 1: export the retained sampling-map image.
# The original map-drawing code and plotted coordinate table are unavailable.

source_map <- here::here("data", "reference", "figures", "Fig01_map.png")
target <- here::here("figures", "main", "Fig01_map.png")
if (!file.exists(source_map)) stop("Missing retained map: ", source_map)
dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
if (!file.copy(source_map, target, overwrite = TRUE)) {
  stop("Could not export Figure 1: ", target)
}
cat(target, "\n")
