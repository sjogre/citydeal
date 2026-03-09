setwd("G:/My Drive/CityDeal/Data_analysis_after_review/")
library(vegan)
library(geosphere)
library(ggplot2)
library(reshape2)
library(dplyr)

# =========================
# INPUT DATA
# =========================
comm_geo <- as.data.frame(phy_rare@otu_table)
comm_geo <- as.data.frame(t(comm_geo))
comm_geo <- comm_geo[order(row.names(comm_geo)), ]


# Community table (samples x taxa)
#comm_geo <- read.csv("", row.names = 1, check.names = FALSE)

# Metadata (must match rownames of comm_geo)
meta_geo <- read.csv("16S_metadata_unrarefied.csv", row.names = 1, dec = ".")
meta_geo <- meta_geo[order(row.names(meta_geo)), ]
rownames(meta_geo)
rownames(comm_geo)


# Ensure sample order matches
meta_geo <- meta_geo[match(rownames(comm_geo), rownames(meta_geo)), ]
meta_geo

# setdiff(rownames(meta_geo), rownames(comm_geo))
 #setdiff(rownames(comm_geo), rownames(meta_geo))
# =========================
# STEP 1: BRAY–CURTIS
# =========================

bc_dist <- vegdist(comm_geo, method = "bray")

# =========================
# STEP 2: GEOGRAPHIC DISTANCE (HAVERSINE)
# =========================

coords <- as.matrix(meta_geo[, c("Latitude", "Longtitude")])

geo_dist <- distm(coords, fun = distHaversine) / 1000
geo_dist <- as.dist(geo_dist)

# =========================
# STEP 3: SOIL TYPE DISTANCE (CATEGORICAL)
# =========================
# 0 = same soil type
# 1 = different soil type

soil <- meta_geo$topology
soil_dist <- as.dist(outer(soil, soil, FUN = "!=") * 1)

# =========================
# STEP 4: PAIRWISE DATAFRAME
# =========================

bc_mat  <- as.matrix(bc_dist)
geo_mat <- as.matrix(geo_dist)

bc_long  <- melt(bc_mat)
geo_long <- melt(geo_mat)

df_geo <- bind_cols(bc_long, geo_long)

colnames(df_geo) <- c("Sample1", "Sample2", "Braycurtis",
                      "Sample3", "Sample4", "GeoDist_km")

df_geo <- df_geo[, c("Sample1", "Sample2", "Braycurtis", "GeoDist_km")]

# Remove self-comparisons and duplicates
head(df_geo)
df_geo <- df_geo[
  df_geo$Sample1 != df_geo$Sample2,
]
head(df_geo)
str(df_geo)
df_geo <- df_geo[as.character(df_geo$Sample1) < as.character(df_geo$Sample2), ]
str(df_geo)
head(df_geo)
# Add soil type info
df_geo$Soil1 <- soil[df_geo$Sample1]
df_geo$Soil2 <- soil[df_geo$Sample2]

# =========================
# STEP 5: DISTANCE–DECAY PLOT (ALL SAMPLES)
# =========================

ggplot(df_geo, aes(x = GeoDist_km, y = Braycurtis)) +
  geom_point(colour = "#4f86f7", alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +
  labs(
    title = "Distance–Decay Relationship (All Soils)",
    x = "Geographic Distance (km)",
    y = "Bray–Curtis Dissimilarity"
  ) +
  theme_bw()

# =========================
# STEP 6: LINEAR MODEL (DESCRIPTIVE ONLY)
# =========================

lm_model <- lm(Braycurtis ~ GeoDist_km, data = df_geo)
summary(lm_model)

# =========================
# STEP 7: PARTIAL MANTEL TEST (KEY CORRECTION)
# =========================
# Tests geography while controlling for soil type

partial_mantel <- mantel.partial(
  bc_dist,
  geo_dist,
  soil_dist,
  method = "pearson",
  permutations = 9999
)
cor(df_geo$GeoDist_km, df_geo$Braycurtis, method = "pearson")
partial_mantel

# =========================
# OPTIONAL: WITHIN-SOIL DISTANCE–DECAY
# =========================
df_within_soil <- df_geo[df_geo$Soil1 == df_geo$Soil2, ]


soil_levels <- unique(df_within_soil$Soil1)


soil_levels <- unique(df_within_soil$Soil1)

# Map sample -> topology using distance matrix rownames
sample_to_topology <- setNames(
  meta_geo$topology,
  meta_geo$SampleID
)

# Build geographic distance matrix
geo_mat <- distm(coords, fun = distHaversine) / 1000

# Attach sample IDs
meta_geo$SampleID <- rownames(meta_geo)
rownames(geo_mat) <- meta_geo$SampleID
colnames(geo_mat) <- meta_geo$SampleID

# Convert to dist object
geo_dist <- as.dist(geo_mat)


library(dplyr)

soil_levels <- unique(df_within_soil$Soil1)

within_soil_stats <- lapply(soil_levels, function(s) {
  print(s)
  # ---- samples that are BOTH in the distance matrix AND this soil
  all_samples <- labels(bc_dist)
  
  samples_s <- all_samples[
    sample_to_topology[all_samples] == s
  ]
  
  samples_s <- samples_s[!is.na(samples_s)]
  n_samples <- length(samples_s)
  print(n_samples)
  # ---- pairwise rows (descriptive stats)
  df_s <- df_within_soil %>% filter(Soil1 == s)
  str(df_s)
  n_pairs <- nrow(df_s)
  print(n_pairs)
  out <- data.frame(
    SoilType = s,
    n_samples = n_samples,
    n_pairs = n_pairs,
    Pearson_r = NA,
    Slope = NA,
    Mantel_r = NA,
    Mantel_p = NA
  )
  
  # ---- descriptive stats (>= 2 pairs)
  if (n_pairs >= 2) {
    out$Pearson_r <- cor(df_s$GeoDist_km, df_s$Braycurtis)
    out$Slope <- coef(lm(Braycurtis ~ GeoDist_km, data = df_s))[2]
  }
  
  # ---- Mantel test (>= 3 samples)
  if (n_samples >= 3) {
    bc_s  <- as.dist(as.matrix(bc_dist)[samples_s, samples_s])
    geo_s <- as.dist(as.matrix(geo_dist)[samples_s, samples_s])
    
    m <- mantel(bc_s, geo_s, permutations = 9999)
    
    out$Mantel_r <- unname(m$statistic)
    out$Mantel_p <- m$signif
  }
  
  out
})

within_soil_stats <- bind_rows(within_soil_stats)
within_soil_stats
within_soil_stats <- bind_rows(within_soil_stats)
within_soil_stats


library(ggplot2)

ggplot(df_within_soil,
       aes(x = GeoDist_km, y = Braycurtis, colour = Soil1)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_manual(
    values = c('#FFABAB', '#ACE7FF', '#DBFFD6', '#ECD4FF', '#FEECC8')
  ) +
  theme_bw() +
  labs(
    title = "Distance–Decay Within Soil Types",
    x = "Geographic Distance (km)",
    y = "Bray–Curtis Dissimilarity",
    colour = "Typology"
  )
ggsave("./ITS_Distance_decay_model.png",width = 15, height = 10, units = "in", dpi = 600)
#ggsave("./S16_Distance_decay_model.png",width = 15, height = 10, units = "in", dpi = 600)
