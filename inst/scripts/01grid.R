### GRID FOR PREDICTION ###

library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("tidyverse")

# projection
ccb_crs <- paste("+proj=aea +lat_1=40",
                 "+lat_2=45 +lat_0=42 +lon_0=-70",
                 "+x_0=0 +y_0=0",
                 "+ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# coords 42.101667, -70.660833 / 42.062333, -70.243028
coordsLimit <- data.frame(matrix(c(-70.660833, -70.243028, 42.101667, 42.062333), 2, 2))
colnames(coordsLimit) <- c("longitude", "latitude")
coordsLimit <- sf::st_as_sf(coordsLimit, coords = c("longitude", "latitude")) %>% 
  sf::st_set_crs(4326) %>% 
  sf::st_transform(ccb_crs)
coordsLimit <- data.frame(st_coordinates(coordsLimit))


# poly is a polygon object of CCB (not publicly available)
grid <- st_make_grid(as(poly, 'sf'), cellsize = 1000, what = "centers")
grid <- st_intersection(grid, as(poly, 'sf'))
grid <- st_coordinates(grid)
m0 <- lm(Y ~ X, data = coordsLimit)
grid <- data.frame(grid[grid[,2] < predict(m0, newdata = data.frame(grid)),])
plot(grid)
dim(grid)

eeuu <- ne_countries(scale = "large", country = "United States of America", returnclass = "sf")
eeuu <- eeuu %>% 
  sf::st_set_crs(4326) %>% 
  sf::st_transform(ccb_crs)

map.ccb <- ggplot(data = eeuu) + 
  geom_sf(fill= "antiquewhite") + 
  coord_sf(xlim = c(-61771.464, 1228.536), ylim = c(-34392.144, 13607.86), expand = FALSE) + 
  xlab("Longitude") + ylab("Latitude") + 
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6,angle=90),
        axis.title=element_text(size=10,face="bold"))
map.ccb

map.ccb <- map.ccb + 
  ggplot2::geom_point(ggplot2::aes(x = X, y = Y), data = grid, color = "grey22")

ggplot2::ggsave("grid.png", map.ccb, width = 8.27 / 2, height = 11.69 / 4)

save(list = c("grid", "eeuu"), file = "grid.RData")
