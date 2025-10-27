# Marina Natividad Avila
# Redoing world Map for GALA project
rm(list=ls())
`%notin%` <- Negate(`%in%`)

library(ggplot2)
theme_set(theme_bw())
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(data.table)
library(tidyr)
library(dplyr)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

coord <- fread('../../worldcities.csv')
coord$key <- paste(coord$city, coord$admin_name, coord$iso2, sep = '::')

coord <- coord[coord$key %in% c('São Paulo::São Paulo::BR', 'Mexico City::Ciudad de México::MX', 'Manhattan::New York::US', 'Miami::Florida::US', 'Oakland::California::US',
                                'Bogotá::Bogotá::CO', 'Lima::Lima::PE', 'Oakland::California::US', 'San José::San José::CR', 
                                'San Diego::California::US', 'San José::San José::CR', 'Kansas City::Missouri::US')]

pedigree <- openxlsx::read.xlsx('../../Supplementary_Tables.xlsx', sheet = 'Supplementary Table 1')
pedigree <- setDT(pedigree)
pedigree <- pedigree[pedigree$Pop == 'AMR']

unique(pedigree$Cohort)
pedigree <- pedigree[pedigree$Cohort %in% c('HERTZ-PICCIOTTO, DAVIS', 'TASC', 'PASSOS-BUENO, BRAZIL', 
                                            'LATTIG, COLOMBIA', 'MAYO, PERU', 'PERICAK-VANCE, MIAMI', 
                                            'KOLEVZON, MSSM', 'CVCR')]

pedigree <- pedigree[, c('Sample', 'Mother', 'Father', 'Cohort', 'VCF_batch', 'Pop')]

long.ped <- pedigree %>% select(Cohort, Sample, Mother, Father) %>%
  tidyr::pivot_longer(cols = c(Sample, Mother, Father), names_to = 'Role', values_to = 'ID')

summary_table <- long.ped %>% distinct(Cohort, ID) %>% 
  group_by(Cohort) %>% summarise(total_unique_IDs = n())

print(summary_table)
#   Cohort                 total_unique_IDs
#    <chr>                             <int>
# 1 CVCR                                545
# 2 HERTZ-PICCIOTTO, DAVIS              439
# 3 KOLEVZON, MSSM                       99
# 4 LATTIG, COLOMBIA                    354
# 5 MAYO, PERU                           96
# 6 PASSOS-BUENO, BRAZIL                634
# 7 PERICAK-VANCE, MIAMI                314
# 8 TASC                                185

coord$label <- c(
  'São Paulo, Brazil N=634',
  'Mexico City, Mexico',
  'Lima, Peru N=96',
  'Bogota, Colombia N=354',
  'Miami, USA N=314',
  'California, USA (CHARGE) N=439',
  'New York, USA N=99',
  'USA, Europe (TASC) N=185',
  'Central Valley, Costa Rica N=545',
  'California, USA (KP)'
)


ggplot(data = world) +
  geom_sf(fill= 'antiquewhite', ) +
  coord_sf(xlim = c(-150, -30), ylim = c(-60, 60), expand = FALSE, datum = NA) +
  theme_void() +
  # annotation_scale(location = "bl", width_hint = 0.25) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +
  ggrepel::geom_label_repel(data = coord, aes(x = lng, y = lat, label = label),
                            size = 4.7, box.padding = 1, point.padding = 0.6, max.overlaps = Inf) +
  geom_point(data = coord, aes(lng, lat), colour = 'black', fill = '#FFC6C6', size = 4, shape = 21) +
  theme(panel.grid.major = element_line(color = gray(.4), linetype = 'blank', size = 0.25), 
        panel.background = element_rect(fill = 'aliceblue'),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.background = element_rect(colour = "black", fill = NA, size = 1))


#ggsave('../../world_map_2025-06-20.pdf')
