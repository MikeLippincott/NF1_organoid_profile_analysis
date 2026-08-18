# set custom colors for each MOA
custom_MOA_palette <- c(
    'BRD4 inhibitor' = "#93152A",  # Dark red
    'receptor tyrosine kinase inhibitor' = "#BA3924",  # Red
    'tyrosine kinase inhibitor' = "#D08543",  # Orange
    'MEK1/2 inhibitor' = "#A1961A",  # Yellow-green/olive

    'IGF-1R inhibitor' = "#9FC62A",  # Yellow-green
    'mTOR inhibitor' = "#1FAD23",  # Green
    'PI3K inhibitor' = "#32D06A",  # Light green
    'PI3K and HDAC inhibitor' = "#15937C",  # Teal/dark green
    'HDAC inhibitor' = "#24A5BA",  # Light blue/cyan

    'Apoptosis induction' = "#438CD0",  # Medium blue
    'DNA binding' = "#1A24A1",  # Dark blue
    'HSP90 inhibitor' = "#532AC6",  # Blue-purple

    'histamine H1\nreceptor antagonist' = "#AD1FA6",  # Purple/magenta
    'Na+/K+ pump inhibitor' = "#D03294",  # Pink/magenta

    'Control' = "#444444"  # Gray
)
treatment_moa_map <- c(
    'DMSO' = 'Control',
    'Staurosporine' = 'Apoptosis induction',
    'Fimepinostat' = 'PI3K and HDAC inhibitor',
    'Copanlisib' = 'PI3K inhibitor',
    'Everolimus' = 'mTOR inhibitor',
    'Rapamycin' = 'mTOR inhibitor',
    'Sapanisertib' = 'mTOR inhibitor',
    'Vistusertib' = 'mTOR inhibitor',
    'Panobinostat' = 'HDAC inhibitor',
    'ARV-825' = 'BRD4 inhibitor',
    'Imatinib' = 'tyrosine kinase inhibitor',
    'Nilotinib' = 'tyrosine kinase inhibitor',
    'Cabozantinib' = 'receptor tyrosine kinase inhibitor',
    'Linsitinib' = 'IGF-1R inhibitor',
    'Binimetinib' = 'MEK1/2 inhibitor',
    'Mirdametinib' = 'MEK1/2 inhibitor',
    'Trametinib' = 'MEK1/2 inhibitor',
    'Selumetinib' = 'MEK1/2 inhibitor',
    "Onalespib" = "HSP90 inhibitor",
    "Digoxin" = "Na+/K+ pump inhibitor",
    "Ketotifen" = "histamine H1\nreceptor antagonist",
    "Trabectedin" = "DNA binding"
)

custom_MOA_order <- c(
    'Control',
    'BRD4 inhibitor',
    'receptor tyrosine kinase inhibitor',
    'tyrosine kinase inhibitor',
    'MEK1/2 inhibitor',
    'IGF-1R inhibitor',
    'mTOR inhibitor',
    'PI3K inhibitor',
    'PI3K and HDAC inhibitor',
    'HDAC inhibitor',
    'Apoptosis induction',
    'DNA binding',
    'HSP90 inhibitor',
    'histamine H1\nreceptor antagonist',
    'Na+/K+ pump inhibitor'
)

# Set custom colors for each treatment
# Hue = MOA family (similar mechanisms share a hue range); Lightness/Saturation = dose
custom_treatment_palette <- c(
    'DMSO' = "#A6A6A6",           # Control

    'Staurosporine' = "#3F468C",  # Broad-spectrum kinase inhibitor

    'Fimepinostat' = "#3DCCA8",   # HDAC/PI3K inhibitor (dual)
    'Copanlisib' = "#3DCACC",     # PI3K inhibitor
    'Everolimus' = "#3DA4CC",     # mTOR inhibitor
    'Rapamycin' = "#3D7DCC",      # mTOR inhibitor
    'Sapanisertib' = "#3D57CC",   # mTOR inhibitor
    'Vistusertib' = "#493DCC",    # mTOR inhibitor

    'Panobinostat' = "#CC8029",   # HDAC inhibitor
    'ARV-825' = "#CCAB29",        # BRD4 inhibitor

    'Imatinib' = "#6047CC",       # BCR-ABL/KIT inhibitor
    'Nilotinib' = "#8347CC",      # BCR-ABL/KIT inhibitor
    'Cabozantinib' = "#A647CC",   # Multi-kinase inhibitor
    'Linsitinib' = "#CA47CC",     # IGF-1R inhibitor

    'Binimetinib' = "#D92B7F",    # MEK inhibitor
    'Mirdametinib' = "#D92B51",   # MEK inhibitor
    'Trametinib' = "#D9342B",     # MEK inhibitor
    'Selumetinib' = "#D9622B",    # MEK inhibitor

    'Onalespib' = "#6CA642",      # HSP90 inhibitor
    'Digoxin' = "#BF3078",        # Na/K-ATPase inhibitor
    'Ketotifen' = "#238C83",      # Antihistamine
    'Trabectedin' = "#388C5B"     # DNA-binding agent
)

dose_palette <- c(
    "1" = "#CFCFCF",
    "10" = "#4D4D4D"
)
custom_treatment_order <- c(
    'DMSO',
    'Staurosporine',
    'Fimepinostat',
    'Copanlisib',
    'Everolimus',
    'Rapamycin',
    'Sapanisertib',
    'Vistusertib',
    'Imatinib',
    'Nilotinib',
    'Cabozantinib',
    'Linsitinib',
    'Panobinostat',
    'ARV-825',
    'Onalespib',
    'Digoxin',
    'Ketotifen',
    'Trabectedin',
    'Binimetinib',
    'Mirdametinib',
    'Trametinib',
    'Selumetinib'
)
channel_palette = c(
    "DNA" = "#0000AB",
    "AGP" = "#b1001a",
    "Mito" = "#B000B0",
    "ER" = "#00D55B",
    "BF" = "#FFFF00",
    "NoChannel" = "#B09FB0"
)

sc_compartment_palette = c(
    "Cell" = "#B000B0",
    "Cytoplasm" = "#00D55B",
    "Nuclei" = "#0000AB"
)

organoid_compartment_palette = c(
    "Organoid" = "#B09FB0"
)

feature_type_palette = c(
    "AreaSizeShape" = brewer.pal(8, "Paired")[1],
    "Colocalization" = brewer.pal(8, "Paired")[2],
    "Granularity" = brewer.pal(8, "Paired")[3],
    "Intensity" = brewer.pal(8, "Paired")[5],
    "Texture" = brewer.pal(8, "Paired")[8]
)
