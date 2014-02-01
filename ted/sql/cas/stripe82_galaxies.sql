--- Yields
--- Rows        kB      Name
--- 592,041     20,552  AllGalsbt20
---

SELECT
    -- objID,
    fieldID,
    psfMag_r,
    run, rerun, camcol, field,
    ra as Ra, dec as Dec

FROM
    Stripe82.PhotoObjAll

WHERE
    run IN (106, 206)
  AND
    -- -9999 is given to undefined (?) values (ref: ??)
    psfMag_r BETWEEN -9999 AND 20.0
  AND
    -- Type 3 means object is categorised as a galaxy in the table.
    type = 3

    -- STRIPE 82 extent
  AND
    (
        Ra BETWEEN 0.0 AND 60.0
      OR
        Ra BETWEEN 300.0 AND 360.0
    )
  AND
   Dec BETWEEN -1.25 AND 1.25

--- ;

