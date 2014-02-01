--- Yields
--- Rows        kB      Name
--- 592,041     20,552  AllGalsbt20
---

SELECT
    objID,
    -- fieldID,
    psfMag_r,
    -- specobjid,
    run, rerun, camcol, field,
    ra as Ra, dec as Dec

    -- Not available in the PhotoObjAll table
    -- raErr, decErr, raDecCorr

FROM
    PhotoObjAll

WHERE
    run IN (106, 206)
  AND
    -- -9999 is given to undefined (?) values (ref: ??)
    psfMag_r BETWEEN -9999 AND 20.0
  AND
    -- Type 3 means object is categorised
    -- as a galaxy in the table.
    type = 3

    -- This takes more than the allowed time to finish (SQL search only--
    -- could not log on to CasJobs at the time).

  --   -- STRIPE 82 extent
  --   -- ----------------
  -- AND
  --   (
  --       -- ra BETWEEN 0 AND 60
  --       (ra >= 0 AND ra <= 60)
  --       -- (0 <= ra <= 60)
  --     OR
  --       -- ra BETWEEN 300 AND 360
  --       (ra >= 300 AND ra <= 360)
  --   )
  -- AND
  --  dec BETWEEN -1.25 AND 1.25
