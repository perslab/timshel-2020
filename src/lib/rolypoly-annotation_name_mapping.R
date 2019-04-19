
# ======================= ***INFO*** ======================= #
# Each data frame MUST set the columns 'name_r', 'name_clean' and 'category'

# ==================================================== #
# ======================= DATA FUNCTION ======================= #
# ==================================================== #

utils.rolypoly_get_data <- function(dataset) {
  if (dataset=="mousebrain") {
    # ======================= Mousebrain ======================= #
    file.data <- here("/data/expression/mousebrain/mousebrain-agg_L5.metadata.csv") # PT formatted/cleaned meta-data
    ### Set new names
    df.category <- suppressMessages(read_csv(file.data)) %>% 
      mutate(
        name_r = annotation,
        name_clean = annotation, # paste0(Class, "-", ClusterName),
        category = Class
      ) 
  } else if (dataset=="gtex") {
    # ======================= GTEX ======================= #
    file.data <- "SMTS	SMTSD	sub_tissue	tissue	category
                                 Adipose Tissue	Adipose - Subcutaneous	Adipose Subcutaneous	Adipose Tissue	Adipose
                                 Adipose Tissue	Adipose - Visceral (Omentum)	Adipose Visceral (Omentum)	Adipose Tissue	Adipose
                                 Adrenal Gland	Adrenal Gland	Adrenal Gland	Adrenal Gland	Other
                                 Blood Vessel	Artery - Aorta	Artery Aorta	Blood Vessel	Cardiovascular
                                 Blood Vessel	Artery - Coronary	Artery Coronary	Blood Vessel	Cardiovascular
                                 Blood Vessel	Artery - Tibial	Artery Tibial	Blood Vessel	Cardiovascular
                                 Bladder	Bladder	Bladder	Bladder	Other
                                 Brain	Brain - Amygdala	Brain Amygdala	Brain	CNS
                                 Brain	Brain - Anterior cingulate cortex (BA24)	Brain Anterior cingulate cortex (BA24)	Brain	CNS
                                 Brain	Brain - Caudate (basal ganglia)	Brain Caudate (basal ganglia)	Brain	CNS
                                 Brain	Brain - Cerebellar Hemisphere	Brain Cerebellar Hemisphere	Brain	CNS
                                 Brain	Brain - Cerebellum	Brain Cerebellum	Brain	CNS
                                 Brain	Brain - Cortex	Brain Cortex	Brain	CNS
                                 Brain	Brain - Frontal Cortex (BA9)	Brain Frontal Cortex (BA9)	Brain	CNS
                                 Brain	Brain - Hippocampus	Brain Hippocampus	Brain	CNS
                                 Brain	Brain - Hypothalamus	Brain Hypothalamus	Brain	CNS
                                 Brain	Brain - Nucleus accumbens (basal ganglia)	Brain Nucleus accumbens (basal ganglia)	Brain	CNS
                                 Brain	Brain - Putamen (basal ganglia)	Brain Putamen (basal ganglia)	Brain	CNS
                                 Brain	Brain - Spinal cord (cervical c-1)	Brain Spinal cord (cervical c-1)	Brain	CNS
                                 Brain	Brain - Substantia nigra	Brain Substantia nigra	Brain	CNS
                                 Breast	Breast - Mammary Tissue	Breast Mammary Tissue	Breast	Other
                                 Blood	Cells - EBV-transformed lymphocytes	Cells EBV-transformed lymphocytes	Blood	Blood/Immune
                                 Skin	Cells - Transformed fibroblasts	Cells Transformed fibroblasts	Skin	Other
                                 Cervix Uteri	Cervix - Ectocervix	Cervix Ectocervix	Cervix Uteri	Other
                                 Cervix Uteri	Cervix - Endocervix	Cervix Endocervix	Cervix Uteri	Other
                                 Colon	Colon - Sigmoid	Colon Sigmoid	Colon	Digestive
                                 Colon	Colon - Transverse	Colon Transverse	Colon	Digestive
                                 Esophagus	Esophagus - Gastroesophageal Junction	Esophagus Gastroesophageal Junction	Esophagus	Digestive
                                 Esophagus	Esophagus - Mucosa	Esophagus Mucosa	Esophagus	Digestive
                                 Esophagus	Esophagus - Muscularis	Esophagus Muscularis	Esophagus	Digestive
                                 Fallopian Tube	Fallopian Tube	Fallopian Tube	Fallopian Tube	Other
                                 Heart	Heart - Atrial Appendage	Heart Atrial Appendage	Heart	Cardiovascular
                                 Heart	Heart - Left Ventricle	Heart Left Ventricle	Heart	Cardiovascular
                                 Kidney	Kidney - Cortex	Kidney Cortex	Kidney	Other
                                 Liver	Liver	Liver	Liver	Liver
                                 Lung	Lung	Lung	Lung	Other
                                 Salivary Gland	Minor Salivary Gland	Minor Salivary Gland	Salivary Gland	Other
                                 Muscle	Muscle - Skeletal	Muscle Skeletal	Muscle	Musculoskeletal/connective
                                 Nerve	Nerve - Tibial	Nerve Tibial	Nerve	Other
                                 Ovary	Ovary	Ovary	Ovary	Other
                                 Pancreas	Pancreas	Pancreas	Pancreas	Pancreas
                                 Pituitary	Pituitary	Pituitary	Pituitary	Other
                                 Prostate	Prostate	Prostate	Prostate	Other
                                 Skin	Skin - Not Sun Exposed (Suprapubic)	Skin Not Sun Exposed (Suprapubic)	Skin	Other
                                 Skin	Skin - Sun Exposed (Lower leg)	Skin Sun Exposed (Lower leg)	Skin	Other
                                 Small Intestine	Small Intestine - Terminal Ileum	Small Intestine Terminal Ileum	Small Intestine	Digestive
                                 Spleen	Spleen	Spleen	Spleen	Blood/Immune
                                 Stomach	Stomach	Stomach	Stomach	Digestive
                                 Testis	Testis	Testis	Testis	Other
                                 Thyroid	Thyroid	Thyroid	Thyroid	Endocrine
                                 Uterus	Uterus	Uterus	Uterus	Other
                                 Vagina	Vagina	Vagina	Vagina	Other
                                 Blood	Whole Blood	Whole Blood	Blood	Blood/Immune"
    ### Set new names
    df.category <- suppressMessages(read_tsv(file.data)) %>% 
      mutate(
        name_r = make.names(SMTSD), # set name_r
        name_clean = SMTSD
      ) 
  } else if (dataset=="maca") {
  # ======================= MACA (v. 180126 - first release) ======================= #
  file.data <- "tissue	cell_type	category
                        Aorta	fibroblast	Cardiovascular
                        Aorta	unknown	Cardiovascular
                        Aorta	epicardial adipocyte	Cardiovascular
                        Aorta	smooth muscle cell	Cardiovascular
                        Aorta	endothelial cell	Cardiovascular
                        Aorta	hematopoietic cell	Cardiovascular
                        Bladder	mesenchymal cell	Other
                        Bladder	bladder cell	Other
                        Bladder	basal cell of urothelium	Other
                        Brain_Microglia	microglial cell	CNS
                        Brain_Microglia	macrophage	CNS
                        Brain_Non-microglia	Bergmann glial cell	CNS
                        Brain_Non-microglia	oligodendrocyte	CNS
                        Brain_Non-microglia	neuron	CNS
                        Brain_Non-microglia	unknown	CNS
                        Brain_Non-microglia	astrocyte of the cerebral cortex	CNS
                        Brain_Non-microglia	oligodendrocyte precursor cell	CNS
                        Brain_Non-microglia	endothelial cell	CNS
                        Brain_Non-microglia	brain pericyte	CNS
                        Brain_Non-microglia	smooth muscle cell	CNS
                        Brain_Non-microglia	neuronal stem cell	CNS
                        Colon	epithelial cell of large intestine	Digestive
                        Colon	large intestine goblet cell	Digestive
                        Colon	enterocyte of epithelium of large intestine	Digestive
                        Colon	enteroendocrine cell	Digestive
                        Colon	Brush cell of epithelium proper of large intestine	Digestive
                        Diaphragm	skeletal muscle satellite stem cell	Other
                        Diaphragm	mesenchymal stem cell	Other
                        Diaphragm	endothelial cell	Other
                        Diaphragm	B cell	Other
                        Diaphragm	macrophage	Other
                        Fat	myeloid cell	Adipose
                        Fat	T cell	Adipose
                        Fat	B cell	Adipose
                        Fat	granulocyte	Adipose
                        Fat	mesenchymal stem cell of adipose	Adipose
                        Fat	endothelial cell	Adipose
                        Fat	natural killer cell	Adipose
                        Fat	epithelial cell	Adipose
                        Fat	neutrophil	Adipose
                        Fat	smooth muscle cell	Adipose
                        Heart	fibroblast	Cardiovascular
                        Heart	endothelial cell	Cardiovascular
                        Heart	smooth muscle cell	Cardiovascular
                        Heart	endocardial cell	Cardiovascular
                        Heart	leukocyte	Cardiovascular
                        Heart	cardiac muscle cell	Cardiovascular
                        Heart	erythrocyte	Cardiovascular
                        Kidney	kidney collecting duct cell	Other
                        Kidney	kidney tubule cell	Other
                        Kidney	endothelial cell	Other
                        Kidney	fenestrated cell	Other
                        Kidney	leukocyte	Other
                        Kidney	fibroblast	Other
                        Liver	natural killer cell	Liver
                        Liver	hepatocyte	Liver
                        Liver	Kupffer cell	Liver
                        Liver	endothelial cell of hepatic sinusoid	Liver
                        Liver	B cell	Liver
                        Lung	endothelial cell	Other
                        Lung	stromal cell	Other
                        Lung	type II pneumocyte	Other
                        Lung	B cell	Other
                        Lung	leukocyte	Other
                        Lung	monocyte	Other
                        Lung	dendritic cell	Other
                        Lung	T cell	Other
                        Lung	Clara cell	Other
                        Lung	ciliated cell	Other
                        Lung	natural killer cell	Other
                        Lung	macrophage	Other
                        Lung	mesothelial cell	Other
                        Lung	epithelial cell	Other
                        Lung	lung neuroendocrine cell	Other
                        Lung	type I pneumocyte	Other
                        Mammary	luminal epithelial cell of mammary gland	Other
                        Mammary	basal cell	Other
                        Mammary	stromal cell	Other
                        Mammary	endothelial cell	Other
                        Marrow	B cell	Blood/Immune
                        Marrow	Fraction A pre-pro B cell	Blood/Immune
                        Marrow	granulocyte	Blood/Immune
                        Marrow	hematopoietic stem cell	Blood/Immune
                        Marrow	monocyte	Blood/Immune
                        Marrow	T cell	Blood/Immune
                        Marrow	natural killer cell	Blood/Immune
                        Marrow	neutrophil	Blood/Immune
                        Muscle	mesenchymal stem cell	Musculoskeletal/connective
                        Muscle	skeletal muscle satellite cell	Musculoskeletal/connective
                        Muscle	endothelial cell	Musculoskeletal/connective
                        Muscle	B cell	Musculoskeletal/connective
                        Muscle	macrophage	Musculoskeletal/connective
                        Muscle	T cell	Musculoskeletal/connective
                        Pancreas	leukocyte	Pancreas
                        Pancreas	pancreatic acinar cell	Pancreas
                        Pancreas	pancreatic PP cell	Pancreas
                        Pancreas	pancreatic ductal cell	Pancreas
                        Pancreas	type B pancreatic cell	Pancreas
                        Pancreas	endothelial cell	Pancreas
                        Pancreas	pancreatic stellate cell	Pancreas
                        Pancreas	pancreatic D cell	Pancreas
                        Pancreas	pancreatic A cell	Pancreas
                        Skin	keratinocyte stem cell	Other
                        Skin	basal cell of epidermis	Other
                        Skin	epidermal cell	Other
                        Skin	stem cell of epidermis	Other
                        Spleen	T cell	Blood/Immune
                        Spleen	B cell	Blood/Immune
                        Spleen	myeloid cell	Blood/Immune
                        Thymus	T cell	Blood/Immune
                        Thymus	mesenchymal stem cell	Blood/Immune
                        Tongue	keratinocyte	Other
                        Tongue	basal cell of epidermis	Other
                        Trachea	stromal cell	Other
                        Trachea	leukocyte	Other
                        Trachea	epithelial cell	Other
                        Trachea	endothelial cell	Other"
  
  
  ### Set new names
  df.category <- suppressMessages(read_tsv(file.data)) %>% 
    mutate(
      name_r = make.names(paste0(tissue, "-", cell_type)), # set name_r | *OBS* only valid for 'per-tissue-celltype' version of MACA
      name_clean = paste0(tissue, "-", cell_type)
    )
  
  # ======================= DEPICT ======================= #
  } else if (dataset=="depict") {
    file.data <- "name_original	category
                                 A10.165.114.830.500.750.Subcutaneous.Fat..Abdominal	Adipose
                                 A10.165.114.830.750.Subcutaneous.Fat	Adipose
                                 A11.329.114.Adipocytes	Adipose
                                 A02.835.583.443.800.800.Synovial.Fluid	Blood/Immune
                                 A10.549.400.Lymph.Nodes	Blood/Immune
                                 A10.549.Lymphoid.Tissue	Blood/Immune
                                 A11.118.637.555.567.562.440.Precursor.Cells..B.Lymphoid	Blood/Immune
                                 A11.118.637.555.567.562.B.Lymphocytes	Blood/Immune
                                 A11.118.637.555.567.569.200.700.T.Lymphocytes..Regulatory	Blood/Immune
                                 A11.118.637.555.567.569.T.Lymphocytes	Blood/Immune
                                 A11.118.637.Leukocytes	Blood/Immune
                                 A11.329.372.600.Macrophages..Alveolar	Blood/Immune
                                 A11.443.Erythroid.Cells	Blood/Immune
                                 A11.627.340.360.Granulocyte.Precursor.Cells	Blood/Immune
                                 A11.627.624.249.Monocyte.Macrophage.Precursor.Cells	Blood/Immune
                                 A11.627.635.Myeloid.Progenitor.Cells	Blood/Immune
                                 A11.872.378.590.635.Granulocyte.Macrophage.Progenitor.Cells	Blood/Immune
                                 A11.872.378.590.817.Megakaryocyte.Erythroid.Progenitor.Cells	Blood/Immune
                                 A11.872.378.Hematopoietic.Stem.Cells	Blood/Immune
                                 A15.145.229.188.Blood.Platelets	Blood/Immune
                                 A15.145.229.637.555.567.562.725.Plasma.Cells	Blood/Immune
                                 A15.145.229.637.555.567.569.200.CD4.Positive.T.Lymphocytes	Blood/Immune
                                 A15.145.229.637.555.Leukocytes..Mononuclear	Blood/Immune
                                 A15.145.229.Blood.Cells	Blood/Immune
                                 A15.145.300.Fetal.Blood	Blood/Immune
                                 A15.145.846.Serum	Blood/Immune
                                 A15.145.Blood	Blood/Immune
                                 A15.378.316.580.Monocytes	Blood/Immune
                                 A15.378.316.Bone.Marrow.Cells	Blood/Immune
                                 A15.382.490.315.583.Neutrophils	Blood/Immune
                                 A15.382.490.555.567.537.Killer.Cells..Natural	Blood/Immune
                                 A15.382.490.555.567.622.Lymphocytes..Null	Blood/Immune
                                 A15.382.490.555.567.Lymphocytes	Blood/Immune
                                 A15.382.520.604.700.Spleen	Blood/Immune
                                 A15.382.520.604.800.Palatine.Tonsil	Blood/Immune
                                 A15.382.680.Phagocytes	Blood/Immune
                                 A15.382.812.260.Dendritic.Cells	Blood/Immune
                                 A15.382.812.522.Macrophages	Blood/Immune
                                 A15.382.812.Mononuclear.Phagocyte.System	Blood/Immune
                                 A15.382.Immune.System	Blood/Immune
                                 A07.231.114.Arteries	Cardiovascular
                                 A07.231.908.670.874.Umbilical.Veins	Cardiovascular
                                 A07.231.908.Veins	Cardiovascular
                                 A07.231.Blood.Vessels	Cardiovascular
                                 A07.541.358.100.Atrial.Appendage	Cardiovascular
                                 A07.541.358.Heart.Atria	Cardiovascular
                                 A07.541.510.110.Aortic.Valve	Cardiovascular
                                 A07.541.560.Heart.Ventricles	Cardiovascular
                                 A07.541.Heart	Cardiovascular
                                 A08.186.211.132.810.428.200.Cerebellum	CNS
                                 A08.186.211.132.Brain.Stem	CNS
                                 A08.186.211.464.405.Hippocampus	CNS
                                 A08.186.211.464.710.225.Entorhinal.Cortex	CNS
                                 A08.186.211.464.Limbic.System	CNS
                                 A08.186.211.653.Mesencephalon	CNS
                                 A08.186.211.730.317.357.352.435.Hypothalamo.Hypophyseal.System	CNS
                                 A08.186.211.730.317.357.Hypothalamus	CNS
                                 A08.186.211.730.317.Diencephalon	CNS
                                 A08.186.211.730.885.287.249.487.Corpus.Striatum	CNS
                                 A08.186.211.730.885.287.249.Basal.Ganglia	CNS
                                 A08.186.211.730.885.287.500.270.Frontal.Lobe	CNS
                                 A08.186.211.730.885.287.500.571.735.Visual.Cortex	CNS
                                 A08.186.211.730.885.287.500.670.Parietal.Lobe	CNS
                                 A08.186.211.730.885.287.500.Cerebral.Cortex	CNS
                                 A08.186.211.865.428.Metencephalon	CNS
                                 A08.186.211.Brain	CNS
                                 A09.371.729.Retina	CNS
                                 A11.872.653.Neural.Stem.Cells	CNS
                                 A03.556.124.369.Intestinal.Mucosa	Digestive
                                 A03.556.124.526.767.Rectum	Digestive
                                 A03.556.124.684.Intestine..Small	Digestive
                                 A03.556.124.Intestines	Digestive
                                 A03.556.249.124.Ileum	Digestive
                                 A03.556.249.249.209.Cecum	Digestive
                                 A03.556.249.249.356.668.Colon..Sigmoid	Digestive
                                 A03.556.249.249.356.Colon	Digestive
                                 A03.556.500.760.464.Parotid.Gland	Digestive
                                 A03.556.500.760.Salivary.Glands	Digestive
                                 A03.556.875.500.Esophagus	Digestive
                                 A03.556.875.875.Stomach	Digestive
                                 A03.556.875.Upper.Gastrointestinal.Tract	Digestive
                                 A03.556.Gastrointestinal.Tract	Digestive
                                 A06.407.071.140.Adrenal.Cortex	Other
                                 A06.407.071.Adrenal.Glands	Other
                                 A06.407.312.782.Testis	Other
                                 A06.407.312.Gonads	Other
                                 A06.407.900.Thyroid.Gland	Other
                                 A06.407.Endocrine.Glands	Other
                                 A10.336.707.Prostate	Other
                                 A11.382.Endocrine.Cells	Other
                                 A11.436.329.Granulosa.Cells	Other
                                 A03.620.Liver	Liver
                                 A11.436.348.Hepatocytes	Liver
                                 A02.165.Cartilage	Musculoskeletal/Connective
                                 A02.633.567.850.Quadriceps.Muscle	Musculoskeletal/Connective
                                 A02.835.232.834.151.Cervical.Vertebrae	Musculoskeletal/Connective
                                 A02.835.583.443.800.Synovial.Membrane	Musculoskeletal/Connective
                                 A05.360.319.679.690.Myometrium	Musculoskeletal/Connective
                                 A10.165.450.300.425.Keloid	Musculoskeletal/Connective
                                 A10.165.450.300.Cicatrix	Musculoskeletal/Connective
                                 A10.690.467.Muscle..Smooth	Musculoskeletal/Connective
                                 A10.690.Muscles	Musculoskeletal/Connective
                                 A11.329.171.Chondrocytes	Musculoskeletal/Connective
                                 A11.329.228.Fibroblasts	Musculoskeletal/Connective
                                 A11.329.629.Osteoblasts	Musculoskeletal/Connective
                                 A11.329.830.Stromal.Cells	Musculoskeletal/Connective
                                 A11.329.Connective.Tissue.Cells	Musculoskeletal/Connective
                                 A11.620.520.Myocytes..Smooth.Muscle	Musculoskeletal/Connective
                                 A03.734.414.Islets.of.Langerhans	Pancreas
                                 A03.734.Pancreas	Pancreas
                                 A04.411.Lung	Other
                                 A04.531.520.Nasal.Mucosa	Other
                                 A05.360.319.114.373.Fallopian.Tubes	Other
                                 A05.360.319.114.630.Ovary	Other
                                 A05.360.319.679.256.Cervix.Uteri	Other
                                 A05.360.319.679.490.Endometrium	Other
                                 A05.360.319.679.Uterus	Other
                                 A05.360.319.887.Vulva	Other
                                 A05.360.319.Genitalia..Female	Other
                                 A05.360.444.492.362.Foreskin	Other
                                 A05.360.444.Genitalia..Male	Other
                                 A05.360.490.Germ.Cells	Other
                                 A05.360.Genitalia	Other
                                 A05.810.453.324.Kidney.Cortex	Other
                                 A05.810.453.Kidney	Other
                                 A05.810.890.Urinary.Bladder	Other
                                 A09.371.Eye	Other
                                 A10.272.497.Epidermis	Other
                                 A10.272.Epithelium	Other
                                 A10.615.284.473.Chorion	Other
                                 A10.615.550.599.Mouth.Mucosa	Other
                                 A10.615.550.Mucous.Membrane	Other
                                 A10.615.789.Serous.Membrane	Other
                                 A10.615.Membranes	Other
                                 A11.436.275.Endothelial.Cells	Other
                                 A11.436.294.064.Glucagon.Secreting.Cells	Other
                                 A11.436.397.Keratinocytes	Other
                                 A11.436.Epithelial.Cells	Other
                                 A11.497.497.600.Oocytes	Other
                                 A11.872.040.Adult.Stem.Cells	Other
                                 A11.872.190.260.Embryoid.Bodies	Other
                                 A11.872.190.Embryonic.Stem.Cells	Other
                                 A11.872.580.Mesenchymal.Stem.Cells	Other
                                 A11.872.700.500.Induced.Pluripotent.Stem.Cells	Other
                                 A11.872.Stem.Cells	Other
                                 A14.549.167.646.Periodontium	Other
                                 A14.549.167.Dentition	Other
                                 A14.549.885.Tongue	Other
                                 A14.549.Mouth	Other
                                 A14.724.557.Nasopharynx	Other
                                 A14.724.Pharynx	Other
                                 A17.815.Skin	Other"
  
    ### DEPICT
    df.category <- read_tsv(file.data)
    ### Cleaning names
    tmp_clean1 <- stringr::str_replace_all(df.category$name_original, pattern="A(.*)(\\d+)\\.", replacement="") # e.g. "Precursor.Cells..B.Lymphoid", "B.Lymphocytes"...
    # tmp_clean1
    tmp_clean2 <- stringr::str_replace_all(tmp_clean1, pattern="\\.\\.", replacement=" - ") # replace double dots ('..') with spaced hypthen (' - ')
    # tmp_clean2
    tmp_clean3 <- stringr::str_replace_all(tmp_clean2, pattern="\\.", replacement=" ") # replace single dots ('.') with space (' ')
    # sort(tmp_clean3)
    ### Set new names
    df.category <- df.category %>% 
      mutate(
        name_r = make.names(paste0(category, "-", name_original)), # set name_r
        name_clean = tmp_clean3
      ) 
  } else {
    stop("Got wrong dataset argument")
  }

  return(df.category)
}





# ======================= Tabula Muris (v. 180920 - final release) ======================= #
# 
# df.category.tabula_muris <- read_tsv("tissue	cell_type	category
#                                      .... <NOT COMPLETED>")
# 
# ### Set new names
# df.category.maca <- df.category.maca %>% 
#   mutate(
#     name_r = make.names(paste0(tissue, "-", cell_type)), # set name_r | *OBS* only valid for 'per-tissue-celltype' version of MACA
#     name_clean = paste0(tissue, "-", cell_type)
#   ) 


