# Cargar librerías ---------------------------------------------------------
require(pacman)
pacman::p_load(terra, sf, fs, tidyverse, rgeos, gtools, stringr, glue, geodata, rnaturalearthdata, rnaturalearth)

# Limpiar memoria y opciones
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999, warn = -1)

# Cargar vector de países (world vector data) ------------------------------
wrld <- ne_countries(returnclass = 'sf', scale = 50)  # Cargar datos vectoriales de países
wrld <- vect(wrld)  # Convertir a formato terra

# Definir ruta local a los archivos descargados ----------------------------
carpeta_datos <- "C:/Users/kuris/OneDrive/Escritorio/tif/co2/TOTALS_emi_nc/"

# Crear función para procesar los archivos ya descargados -------------------
procesar_archivos <- function(yr){
  
  ### Mensaje para indicar el año que se procesa
  cat('Procesando el año ', yr, '\n')
  
  ### Ruta del archivo .nc correspondiente al año
  archivo_nc <- glue("{carpeta_datos}{yr}_TOTALS_emi.nc")
  
  ### Verificar si el archivo existe
  if (file.exists(archivo_nc)) {
    # Aquí puedes procesar el archivo como desees
    cat('¡Archivo encontrado y listo para procesar!\n')
  } else {
    cat('¡Archivo no encontrado!\n')
  }
}

## Aplicar la función a los años deseados ----------------------------------
map(1970:2022, procesar_archivos)

# Crear un raster único y extraer por máscara para el mundo ----------------
fles <- dir_ls(carpeta_datos, regexp = '.nc$')  # Listar los archivos .nc
rstr <- rast(fles)  # Leer todos los archivos como un solo objeto raster
rstr <- terra::crop(rstr, wrld)  # Recortar el raster al área del mundo
rstr <- terra::mask(rstr, wrld)  # Aplicar la máscara de los países

# Extraer el año de los nombres de archivos y asignar nombres a las capas ---
year <- basename(fles) %>% str_split(., pattern = '_') %>% map_chr(5)
names(rstr) <- glue('emi_ton_{year}')  # Asignar nombres de capa por año

# Definir la ruta para guardar el archivo
output_path <- "C:/Users/kuris/OneDrive/Escritorio/tif/co2/tif/"

# Crear el directorio si no existe
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Guardar el raster de la serie temporal -----------------------------------
terra::writeRaster(x = rstr, filename = glue('{output_path}emi_tons_timeseries.tif'), overwrite = TRUE)

# Functions ---------------------------------------------------------------
theme_for_the_win <- function(){
  theme_void() +
    theme(
      legend.position = "top",
      legend.title = element_text(
        size = 9, color = "grey20"
      ),
      legend.text = element_text(
        size = 9, color = "grey20"
      ),
      plot.margin = unit(
        c(
          t = 1, r = 0, # Add 1
          b = 0, l = 0 
        ), "lines"
      )
    )
}

# World vector data -------------------------------------------------------
wrld <- ne_countries(returnclass = 'sf', scale = 50)

# Raster data -------------------------------------------------------------
rstr <- terra::rast('C:/Users/kuris/OneDrive/Escritorio/tif/co2/tif/emi_tons_timeseries.tif')

library(exactextractr)
# To calculate the zonal  -------------------------------------------------
znal <- exact_extract(x = rstr, y = wrld, fun = 'sum')
znal <- as_tibble(znal)
znal <- mutate(znal, iso = wrld$iso_a3)
znal <- dplyr::select(znal, iso, everything())

# All the values as a vector ----------------------------------------------
vles <- round(as.numeric(unlist(as.vector(znal[,2:ncol(znal)]))), 0)
brks <- classInt::classIntervals(var = vles, n = 5, style = 'fisher')
brks <- brks$brks
brks <- round(brks, -2.5)
clss <- tibble(class = 1:5, min = brks[1:5], max = brks[2:6], interval = glue('{min}-{max}'))

# To join the table with the shapefile ------------------------------------
shpf <- inner_join(wrld, znal, by = c('iso_a3' = 'iso'))

# To make the map  --------------------------------------------------------
yr=2009
## Function to use -----------
mke.map <- function(yr){
  
  cat('>>> Process: ', yr, '\n')
  shp <- dplyr::select(shpf, iso_a3, name, glue('sum.emi_ton_{yr}'), geometry)
  colnames(shp)[3] <- 'value'
  vls <- pull(shp, value)
  
  # To find the interval for the values
  shp <- mutate(shp, class = findInterval(x = vls, vec = brks, all.inside = TRUE))
  shp <- inner_join(shp, clss, by = 'class')
  shp <- mutate(shp, interval = factor(interval, levels = clss$interval))
  
  # To make the map
  gmp <- ggplot() + 
    geom_sf(data = shp, aes(fill = interval)) + 
    scale_fill_manual(values = brewer.pal(n = 5, name = 'YlOrRd')) + 
    coord_sf() + 
    theme_minimal() + 
    labs(fill = 'GEI (Ton)', 
         caption = 'Adaptado de EDGAR (Emissions Database for Global Atmospheric Research)') + 
    ggtitle(label = glue('Emisiones totales de GEI por país - {yr}')) +
    theme(
      legend.position = 'bottom', 
      text = element_text(family = 'Segoe UI'), 
      plot.title = element_text(size = 16, face = 'bold', hjust = 0.5)
    ) +
    guides(fill = guide_legend( 
      direction = 'horizontal',
      keyheight = unit(3.5, units = "mm"),
      keywidth = unit(35, units = "mm"),
      title.position = 'top',
      title.hjust = 0.5,
      label.hjust = .5,
      nrow = 1,
      byrow = T,
      reverse = F,
      label.position = "bottom"
    )) 
  gmp
  
  ggsave(plot = gmp, filename = glue('./png/maps/gei_{yr}.jpg'), units = 'in', width = 9, height = 6.5, dpi = 300)
  
  # Finish!
  cat('Done!\n')
  
}

-----------------------------------------------------------------------
print(rstr)

layer_index <- year - 1970 + 1  # Asegúrate de que indexas desde 1
if (layer_index <= nlyr(rstr)) {
  plot(rstr[[layer_index]], col = brewer.pal(9, "YlGnBu"), main = paste("Emisiones de CO2 en", year))
} else {
  cat("No hay datos para el año:", year, "\n")
}

mke.map <- function(year) {
  year <- as.numeric(year)  # Asegurarse de que 'year' es numérico
  cat("Creando mapa para el año:", year, "\n")
  
  # Índice de la capa correspondiente al año
  layer_index <- year - 1970 + 1
  
  # Verificar si el índice está dentro del rango de capas en el raster
  if (layer_index <= nlyr(rstr)) {
    plot(rstr[[layer_index]], col = brewer.pal(9, "YlGnBu"), main = paste("Emisiones de CO2 en", year))
  } else {
    cat("No hay datos para el año:", year, "\n")
  }
}

# Aplicar la función para los años 1970 a 2022
map(1970:2022, mke.map)


# To make the GIF from the maps -------------------------------------------
# Cargar las librerías necesarias para crear el GIF
library(magick)
library(fs)

# Crear el GIF a partir de los mapas PNG generados
imgs <- dir_ls('./png/maps', regexp = '.png$')  # Listar los archivos PNG
imgs <- map(imgs, image_read)  # Leer los PNG
jnds <- image_join(imgs)  # Unir las imágenes en una secuencia
anmt <- magick::image_animate(jnds, fps = 10)  # Crear la animación (GIF)

# Crear la carpeta 'gif' si no existe y guardar el GIF
dir_create('./gif')
image_write(image = anmt, path = './gif/gei_10fps.gif')  # Guardar el GIF

## To compile the GIF
imgs <- dir_ls('./png/maps')
imgs <- map(imgs, image_read)
jnds <- image_join(imgs)
anmt <- magick::image_animate(jnds, fps = 10)

dir_create('./gif')

## To write the GIF
image_write(image = anmt, path = './gif/gei_10fps.gif')

---------------------------------------------------

 