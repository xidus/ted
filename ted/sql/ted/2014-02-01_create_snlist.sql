/*
 * SQLite CREATE TABLE sne
 * table for the Supernovae (SNe)
 *
 */

CREATE TABLE Supernovae (
    id          INTEGER     PRIMARY KEY ASC,
    SDSS_id     VARCHAR(7)  NOT NULL,
    SN_type     VARCHAR(15) ,
    IAUC_id     VARCHAR(10) ,
    Ra          REAL        NOT NULL,
    Dec         REAL        NOT NULL,
    redshift    REAL        ,
    Peak_MJD    REAL        --,
--    Flag        VARCHAR(5)  NOT NULL DEFAULT ''
);

--    run     SMALLINT    NOT NULL,
--    rerun   SMALLINT    NOT NULL,
--    camcol  SMALLINT    NOT NULL,
--    field   SMALLINT    NOT NULL,
--
--    ra      REAL        NOT NULL,
--    dec     REAL        NOT NULL

--    FITS_HEADER_01 CHAR(20) NOT NULL,
--    FITS_HEADER_02 CHAR(20) NOT NULL,
--    FITS_HEADER_03 CHAR(20) NOT NULL,
--    FITS_HEADER_04 CHAR(20) NOT NULL,
--    FITS_HEADER_05 CHAR(20) NOT NULL


