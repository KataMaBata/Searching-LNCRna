install.packages("writexl")

library(rvest)
library(httr)
library(xml2)
library(dplyr)
library(stringr)
library(writexl)

setwd("C:/Users/yerbi/Desktop/turbinator przeszukiwator v1.0")

# Główny adres FTP
base_url <- "http://ftp.ensemblgenomes.org/pub/fungi/release-61/fasta/"

# Pobierz listę organizmów (folderów)
get_dirs <- function(url) {
  page <- read_html(url)
  links <- html_nodes(page, "a") %>% html_attr("href")
  dirs <- links[str_detect(links, "/$") & links != "../"]
  return(dirs)
}

organism_dirs <- get_dirs(base_url)

# Funkcja do pobrania informacji o pliku
get_file_size <- function(file_url) {
  response <- HEAD(file_url)
  if (status_code(response) == 200) {
    size <- as.numeric(headers(response)[["content-length"]])
    return(size)
  }
  return(NA)
}

# Lista wyników
results <- list()

# Przejście przez każdy organizm
for (org in organism_dirs) {
  org_url <- paste0(base_url, org)
  
  # Ścieżki do ncrna i cds
  ncrna_url <- paste0(org_url, "ncrna/")
  cds_url <- paste0(org_url, "cds/")
  
  # Pobierz pliki z katalogu ncrna/
  ncrna_files <- tryCatch(get_dirs(ncrna_url), error = function(e) NULL)
  cds_files <- tryCatch(get_dirs(cds_url), error = function(e) NULL)
  
  # Pobierz nazwę pliku ncrna
  ncrna_page <- tryCatch(read_html(ncrna_url), error = function(e) NULL)
  cds_page <- tryCatch(read_html(cds_url), error = function(e) NULL)
  
  if (!is.null(ncrna_page)) {
    ncrna_links <- html_nodes(ncrna_page, "a") %>% html_attr("href")
    ncrna_file <- ncrna_links[str_detect(ncrna_links, "\\.ncrna\\.fa\\.gz$")]
    if (length(ncrna_file) > 0) {
      ncrna_full_url <- paste0(ncrna_url, ncrna_file[1])
      ncrna_size <- get_file_size(ncrna_full_url)
    } else {
      ncrna_full_url <- NA
      ncrna_size <- NA
    }
  } else {
    ncrna_full_url <- NA
    ncrna_size <- NA
  }
  
  if (!is.null(cds_page)) {
    cds_links <- html_nodes(cds_page, "a") %>% html_attr("href")
    cds_file <- cds_links[str_detect(cds_links, "\\.cds\\.all\\.fa\\.gz$")]
    if (length(cds_file) > 0) {
      cds_full_url <- paste0(cds_url, cds_file[1])
      cds_size <- get_file_size(cds_full_url)
    } else {
      cds_full_url <- NA
      cds_size <- NA
    }
  } else {
    cds_full_url <- NA
    cds_size <- NA
  }
  
  # Dodaj dane do listy
  results[[length(results) + 1]] <- data.frame(
    organism = gsub("/$", "", org),
    ncrna_url = ncrna_full_url,
    ncrna_size = ncrna_size,
    cds_url = cds_full_url,
    cds_size = cds_size,
    total_size = sum(c(ncrna_size, cds_size), na.rm = TRUE)
  )
}


# Połącz dane i odfiltruj organizmy bez ncrna i cds
final_df <- bind_rows(results) %>%
  filter(!is.na(ncrna_url)) %>%  # ← Zachowaj tylko te organizmy, które mają ncRNA
  mutate(
    organism_clean = str_replace(organism, "_\\d+.*$", ""),
    organism_clean = str_replace_all(organism_clean, "_", " "),
    organism_clean = str_to_title(organism_clean)
  ) %>%
  arrange(desc(ncrna_size)) %>%
  select(organism_clean, organism, ncrna_url, ncrna_size, cds_url, cds_size, total_size)



# Wyświetl top 10 największych
print(head(final_df, 10))

write.csv2(final_df, "ensembl_fungi_list", row.names = FALSE)
write_xlsx(final_df, "ensembl_fungi_list.xlsx")


# Pobieranie ncRNA oraz cds dla każdego organizmu-----------

# Katalog roboczy
output_dir <- getwd()

# Funkcja do pobierania pliku
download_file <- function(url, destfile) {
  if (!is.na(url)) {
    tryCatch({
      download.file(url, destfile, mode = "wb")
      message("Pobrano: ", destfile)
    }, error = function(e) {
      message("Błąd pobierania: ", url)
    })
  }
}

# Pobierz pliki dla każdego organizmu
for (i in seq_len(nrow(final_df))) {
  org_row <- final_df[i, ]
  
  # Nazwa folderu na podstawie organism_clean
  org_folder <- file.path(output_dir, org_row$organism_clean)
  dir.create(org_folder, showWarnings = FALSE)
  
  # Pobierz plik ncRNA
  if (!is.na(org_row$ncrna_url)) {
    ncRNA_filename <- paste0(org_row$organism, ".ncrna.fa.gz")
    ncRNA_dest <- file.path(org_folder, ncRNA_filename)
    download_file(org_row$ncrna_url, ncRNA_dest)
  }
  
  # Pobierz plik CDS (jeśli istnieje)
  if (!is.na(org_row$cds_url)) {
    cds_filename <- paste0(org_row$organism, ".cds.all.fa.gz")
    cds_dest <- file.path(org_folder, cds_filename)
    download_file(org_row$cds_url, cds_dest)
  }
}

