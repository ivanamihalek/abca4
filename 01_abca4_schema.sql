
-- `eye` varchar(10) ENUM('better', 'average', 'both') = better of OS, OD, average,
--  or the measurement taken with both eyes open (OU)
--  disease progression is recorded as age1:value1;age2:value2
--  patient id refers to the id in the publication
DROP TABLE IF EXISTS `meta`;
CREATE TABLE `meta` (
  `id` mediumint(9) NOT NULL AUTO_INCREMENT,
  `name` varchar(50),
  `value` text,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `cases`;
CREATE TABLE `cases` (
  `id` mediumint(9) NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `allele_id_1`   mediumint(9),
  `allele_id_2`   mediumint(9),
  `publication_id`    int,
  -- this refers to the id in the original publication
  `patient_xref_id` varchar(50),
  `onset_age` int,
  `acuity_type` varchar(50),
  `eye` ENUM ('better', 'average', 'both'),
  `progression` text
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `alleles`;
CREATE TABLE `alleles` (
  `id` mediumint(9) NOT NULL,
  `variant_id` mediumint(9) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


-- int range -2,147,483,648 to 2,147,483,647
-- chrom 1 in human 248,956,422 bp (ABCA4 is on chrom 1)
-- subsitution is referred for a single nucleotide substitution
-- delins for more than 1 or for inversion
DROP TABLE IF EXISTS `variants`;
CREATE TABLE `variants`(
  `id` mediumint(9) NOT NULL AUTO_INCREMENT PRIMARY KEY,
  `gdna_start` int NOT NULL,
  `gdna_end`   int NOT NULL,
  `mod_type` ENUM ('del', 'ins', 'sub', 'delins'),
  `nt_from`  varchar(50),
  `nt_to`    varchar(50),
  `cdna`    varchar(100),
  `cdna_exp_verification` varchar(100) ,
  `cdna_exp_publication` mediumint(9),
  `protein` varchar(100),
  `protein_exp_verification` varchar(100),
  `protein_exp_publication` mediumint(9),
  -- loop: loop between domains
  `protein_domain` ENUM ('NBD1', 'NBD2', 'TMD1', 'TMD2', 'ECD1', 'ECD2', 'R', 'N-term', 'C-term' ,'loop'),
  -- null refers to no protein product
  -- ? means unknown
  -- [membrane, splicing, folding, transport]: impaired [membrane incorporation, splicing, folding, payload transport]
  `systems_effect` ENUM ('null', '?', 'membrane', 'splicing', 'folding', 'transport') default '?',
  `gnomad_freq` float,
  `gnomad_homozygotes` smallint
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



DROP TABLE IF EXISTS `publications`;
CREATE TABLE `publications` (
   `id` mediumint(9) NOT NULL AUTO_INCREMENT  PRIMARY KEY ,
   `pubmed`   int,
   `reference`  text CHARACTER SET utf8mb4 NOT NULL,
   `pubmedcentral`    varchar(20),
   `other_xref` text
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
