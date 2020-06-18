
-- `eye` varchar(10) ENUM('better', 'average', 'both') = better of OS, OD, average,
--  or the measurement taken with both eyes open (OU)
--  disease progression is recorded as age1:value1;age2:value2
--  patient id refers to the id in the publication
DROP TABLE IF EXISTS `cases`;
CREATE TABLE `cases` (
  `id` mediumint(9) NOT NULL AUTO_INCREMENT,
  `allele1_cdna`      varchar(100),
  `allele1_protein`   varchar(50),
  `allele2_cdna`      varchar(100),
  `allele2_protein`   varchar(50),
  `pubmed`    int,
  `patient_id` varchar(50),
  `acuity_type` varchar(50),
  `eye` ENUM('better', 'average', 'both'),
  `progression` text,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `publications`;
CREATE TABLE `publications` (
   `id` mediumint(9) NOT NULL AUTO_INCREMENT,
   `pubmed`   int,
   `reference`  text CHARACTER SET utf8mb4 NOT NULL,
   `pubmedcentral`    varchar(20),
   `other_xref` text,
    PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
