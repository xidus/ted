/*
 * SQLite CREATE TABLE snlist
 * table for the Supernovae (SNe)
 *
 */

CREATE TABLE snlist (
    id          INTEGER     PRIMARY KEY ASC,
    SDSS_id     VARCHAR(10) NOT NULL,
    SN_type     VARCHAR(15) ,
    IAUC_id     VARCHAR(10) ,
    Ra          REAL        NOT NULL,
    Dec         REAL        NOT NULL,
    redshift    REAL        ,
    Peak_MJD    REAL
);

