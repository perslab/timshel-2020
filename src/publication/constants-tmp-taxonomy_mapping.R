
df.taxonomy_metadata <- read_tsv("TaxonomyRank1	TaxonomyRank2	TaxonomyRank3	TaxonomyRank4	TaxonomyRank4_merged1	order_idx_table_s5	order_idx_mb_taxonomy_fig1c	N_CellTypes	CellTypes
Immune cells			Perivascular macrophages	Immune/Microglia	38	1	2	PVM2,PVM1
Immune cells			Microglia	Immune/Microglia	39	1	3	MGL1,MGL2,MGL3
Vascular cells			Vascular and leptomeningeal cells	Vascular	34	2	3	VLMC1,ABC,VLMC2
Vascular cells			Vascular smooth muscle cells	Vascular	35	2	3	VSMCA,VECA,PER3
Vascular cells			Pericytes	Vascular	36	2	2	PER2,PER1
Vascular cells			Vascular endothelial cells	Vascular	37	2	2	VECV,VECC
Glia	Neural crest-like glia		Enteric glia	Enteric glia	33	3	9	ENTG2,ENTG3,ENTG4,ENTG1,ENTG5,ENTG6,ENTG7,ENMFB
Glia	Neural crest-like glia		Olfactory ensheathing cells	Other glia	29	4	1	OEC
Glia	Neural crest-like glia		Schwann cells	Other glia	31	4	1	SCHW
Glia	Neural crest-like glia		Satellite glia	Other glia	32	4	2	SATG2,SATG1
Glia	CNS glia	Astroependymal cells	Astrocytes	Astrocytes	28	5	7	ACNT2,ACNT1,ACOB,ACMB,ACTE2,ACTE1,ACBG
Glia	CNS glia	Astroependymal cells	Choroid epithelial cells	Other CNS glia	23	6	1	CHOR
Glia	CNS glia	Astroependymal cells	Subcommissural organ hypendymal cells	Other CNS glia	24	6	1	HYPEN
Glia	CNS glia	Astroependymal cells	Ependymal cells	Other CNS glia	25	6	3	EPSC,EPEN,EPMB
Glia	CNS glia	Astroependymal cells	Dentate gyrus radial glia-like cells	Other CNS glia	26	6	1	RGDG
Glia	CNS glia	Astroependymal cells	Subventricular zone radial glia-like cells	Other CNS glia	27	6	1	RGSZ
Glia	CNS glia	Oligodendrocytes	Oligodendrocytes	Oligodendrocytes	22	7	9	NFOL1,NFOL2,MFOL2,MOL2,COP2,MOL1,MFOL1,COP1,MOL3
Glia	Neural crest-like glia		Oligodendrocyte precursor cells	Oligodendrocytes	30	7	1	OPC
Neurons	PNS neurons	Peripheral sensory neurons	Peripheral sensory peptidergic neurons	Peripheral sensory neurons	19	8	8	PSPEP1,PSPEP4,PSPEP2,PSPEP3,PSPEP5,PSPEP8,PSPEP6,PSPEP7
Neurons	PNS neurons	Peripheral sensory neurons	Peripheral sensory neurofilament neurons	Peripheral sensory neurons	20	8	3	PSNF1,PSNF3,PSNF2
Neurons	PNS neurons	Peripheral sensory neurons	Peripheral sensory non-peptidergic neurons	Peripheral sensory neurons	21	8	6	PSNP5,PSNP6,PSNP1,PSNP4,PSNP3,PSNP2
Neurons	PNS neurons	Sympathetic neurons	Sympathetic noradrenergic neurons	Sympathetic neurons	17	9	5	SYNOR2,SYNOR1,SYNOR5,SYNOR3,SYNOR4
Neurons	PNS neurons	Sympathetic neurons	Sympathetic cholinergic neurons	Sympathetic neurons	18	9	2	SYCHO2,SYCHO1
Neurons	PNS neurons	Enteric neurons	Enteric neurons	Enteric neurons	16	10	9	ENT9,ENT8,ENT1,ENT5,ENT6,ENT2,ENT4,ENT7,ENT3
Neurons	CNS neurons	Hindbrain neurons	Hindbrain neurons	Hindbrain neurons	11	12	15	HBINH3,HBGLU9,HBGLU5,HBINH8,HBINH4,HBCHO2,HBINH7,HBINH2,HBINH6,HBGLU7,HBGLU6,HBCHO1,HBGLU8,HBGLU4,SCINH1
Neurons	CNS neurons	Cerebellum neurons	Cerebellum neurons	Hindbrain neurons	12	12	6	CBGRC,CBINH1,MEINH1,CBINH2,CBNBL2,CBPC
Neurons	CNS neurons	Di- and mesencephalon neurons	Di- and mesencephalon inhibitory neurons	Di- and mesencephalon inhibitory neurons	10	13	17	DEINH3,MEINH2,MEINH3,TEINH1,HBINH5,MEINH13,MEINH7,MEINH4,MEINH12,MEINH9,MEINH6,MEINH10,MEINH11,MEINH8,HBINH1,MEINH5
Neurons	CNS neurons	Di- and mesencephalon neurons	Di- and mesencephalon excitatory neurons	Di- and mesencephalon excitatory neurons	9	14	22	MEGLU1,DEGLU5,MEGLU10,MEGLU11,DEGLU4,HBGLU3,MEGLU7,MEGLU8,MEGLU2,MEGLU3,HBGLU2,MEGLU6,MEGLU9,MBCHO1,DEGLU1,MEGLU4,MEGLU5,HBGLU1,DEGLU3,DECHO2,DEGLU2,CR
Neurons	CNS neurons	Spinal cord neurons	Spinal cord excitatory neurons	Spinal cord excitatory neurons	7	15	11	SCGLU3,SCGLU4,SCGLU6,SCGLU2,SCGLU7,SCGLU1,SCGLU10,HBGLU10,SCGLU9,SCGLU5,SCGLU8
Neurons	CNS neurons	Spinal cord neurons	Spinal cord inhibitory neurons	Spinal cord inhibitory neurons	8	16	10	HBINH9,SCINH10,SCINH8,SCINH2,SCINH4,SCINH7,SCINH3,SCINH9,SCINH5,SCINH6
Neurons	CNS neurons	Cholinergic, monoaminergic and peptidergic neurons	Peptidergic neurons	Peptidergic neurons	6	17	15	DEINH5,DEINH6,HYPEP1,DEINH7,HYPEP3,TEINH3,HYPEP5,TEINH2,DEINH4,HYPEP4,HYPEP2,MEINH14,DEINH8,SCINH11,HYPEP8
Neurons	CNS neurons	Cholinergic, monoaminergic and peptidergic neurons	Cholinergic and monoaminergic neurons	Cholinergic and monoaminergic neurons	5	18	16	HBSER5,HBCHO4,HBNOR,TECHO,HBCHO3,HBSER3,DECHO1,MEGLU14,MBDOP1,HBSER4,HBSER2,HBSER1,HBADR,HYPEP7,MBDOP2,HYPEP6
Neurons	CNS neurons	Telencephalon interneurons	Telencephalon inhibitory interneurons	Telencephalon inhibitory interneurons	4	19	20	TEINH12,TEINH11,TEINH17,TEINH15,TEINH18,TEINH9,TEINH4,TEINH13,TEINH19,TEINH5,DEINH2,TEINH16,TEINH10,DEINH1,TEINH7,TEINH8,TEINH20,TEINH14,TEINH21,TEINH6
Neurons	CNS neurons	Telencephalon interneurons	Olfactory inhibitory neurons	Olfactory inhibitory neurons	3	20	9	OBINH5,OBINH4,OBDOP1,OBINH1,OBNBL5,OBINH3,OBNBL4,OBINH2,OBDOP2
Neurons	CNS neurons	Immature neural	Glutamatergic neuroblasts	Neuroblasts	1	21	4	OBNBL2,OBNBL1,SEPNBL,CBNBL1
Neurons	CNS neurons	Immature neural	Non-glutamatergic neuroblasts	Neuroblasts	2	21	5	DGNBL2,SZNBL,OBNBL3,DGNBL1,DETPH
Neurons	CNS neurons	Telencephalon projecting neurons	Dentate gyrus granule neurons	Other CNS neurons	14	22	2	DGGRC2,DGGRC1
Neurons	CNS neurons	Telencephalon projecting neurons	Telencephalon projecting inhibitory neurons	Telencephalon projecting inhibitory neurons	15	23	6	MSN1,MSN5,MSN2,MSN4,MSN3,MSN6
Neurons	CNS neurons	Telencephalon projecting neurons	Telencephalon projecting excitatory neurons	Telencephalon projecting excitatory neurons	13	24	24	TEGLU23,TEGLU17,TEGLU4,TEGLU21,TEGLU24,TEGLU19,TEGLU9,TEGLU3,TEGLU6,TEGLU8,TEGLU11,TEGLU7,TEGLU2,TEGLU10,TEGLU5,TEGLU20,TEGLU15,TEGLU13,TEGLU18,TEGLU12,TEGLU22,TEGLU1,TEGLU16,TEGLU14")
