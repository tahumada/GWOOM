"""
BOOM alert filters used for GW (EM-counterpart) follow-up.

Each filter is a MongoDB aggregation pipeline (a list of stages) that BOOM
runs on every incoming alert. An alert that makes it through all stages
passes the filter. Every pipeline has the same layout:

    1. $match      - cheap cuts on alert-packet fields (real/bogus, artifacts,
                     moving objects, known stars)
    2. $project    - keep only the fields the later stages need
    3. $addFields  - compute boolean flags from catalog cross-matches
                     (stars, AGN/QSOs, nuclear sources)
    4. $match      - reject on those flags (+ age / light-curve cuts)
    5. $project    - choose which fields the output keeps

Each filter is written out in full so it can be read, copied and edited on
its own.

Filters:
    ZTF_BASIC              - young, real ZTF transients; not stars or QSOs
    LSST_EMGW_NON_NUCLEAR  - young LSST transients; not stars or AGN/QSOs
    LSST_EMGW_BASIC        - LSST_EMGW_NON_NUCLEAR + reject sources on a galaxy nucleus
    LSST_EMGW_FAST         - LSST_EMGW_BASIC + fast fading/rising in at least one band
"""


# =============================================================================
# ZTF_BASIC
# Young (< 30 d), real, non-moving ZTF transients that are not known stars
# and have no Milliquas quasar counterpart.
# =============================================================================

ZTF_BASIC = [
    # --- 1. Alert quality and contaminant cuts ---
    {"$match": {"$and": [
        {"candidate.isdiffpos": {"$eq": True}},           # positive difference flux (brightening)
        {"candidate.rb": {"$gte": 0.2}},                  # real/bogus score
        {"candidate.drb": {"$gte": 0.8}},                 # deep-learning real/bogus score
        {"properties.rock": {"$eq": False}},              # not a known asteroid
        {"properties.stationary": {"$eq": True}},         # not moving between detections
        {"properties.star": {"$eq": False}},              # not a known star
        {"properties.near_brightstar": {"$eq": False}},   # not in a bright-star halo
        {"candidate.ndethist": {"$gte": 1}},              # at least one detection in history
    ]}},

    # --- 2. Keep only the fields used below ---
    {"$project": {
        "objectId": 1,
        "candidate.jd": 1,
        "candidate.isdiffpos": 1,
        "candidate.rb": 1,
        "candidate.drb": 1,
        "properties.rock": 1,
        "properties.stationary": 1,
        "properties.star": 1,
        "properties.near_brightstar": 1,
        "candidate.ndethist": 1,
        "candidate.jdstarthist": 1,
        "cross_matches.milliquas_v8": 1,
    }},

    # --- 3. Flag alerts with a Milliquas match that is a quasar ---
    # (Q = QSO, R = radio, X = X-ray, 2 = 2QZ/6QZ)
    {"$addFields": {
        "Is_Quasar": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.milliquas_v8", []]},
            "in": {"$in": ["$$this.objtype", ["Q", "Q2", "Q2X", "QR", "QR2", "QR2X", "QRX", "QX"]]},
        }}},
    }},

    # --- 4. Not a quasar, and first detected at most 30 days ago ---
    {"$match": {"$and": [
        {"Is_Quasar": {"$eq": False}},
        {"$expr": {"$lte": [{"$subtract": ["$candidate.jd", "$candidate.jdstarthist"]}, 30]}},
    ]}},

    # --- 5. Output fields ---
    {"$project": {
        "objectId": 1,
        "candidate.jd": 1,
        "candidate.isdiffpos": 1,
        "candidate.rb": 1,
        "candidate.drb": 1,
        "properties.rock": 1,
        "properties.stationary": 1,
        "properties.star": 1,
        "properties.near_brightstar": 1,
        "candidate.ndethist": 1,
        "candidate.jdstarthist": 1,
        "cross_matches.milliquas_v8": 1,
        "Is_Quasar": 1,
    }},
]


# =============================================================================
# LSST_EMGW_NON_NUCLEAR
# Young (< 30 d), clean, non-moving LSST transients that are not stars
# (Gaia / LSPSC / PS1) and not AGN/QSOs (WISE colour / Milliquas).
# Unlike LSST_EMGW_BASIC, there is no veto on sources at a galaxy nucleus.
# =============================================================================

LSST_EMGW_NON_NUCLEAR = [
    # --- 1. Alert quality and contaminant cuts ---
    {"$match": {"$and": [
        {"candidate.isDipole": {"$eq": False}},           # not a dipole (bad subtraction)
        {"candidate.glint_trail": {"$eq": False}},        # not a satellite glint / trail
        {"candidate.isdiffpos": {"$eq": True}},           # positive difference flux (brightening)
        {"properties.rock": {"$eq": False}},              # not a known asteroid
        {"properties.stationary": {"$eq": True}},         # not moving between detections
        {"properties.star": {"$eq": False}},              # not a known star
        {"properties.near_brightstar": {"$eq": False}},   # not in a bright-star halo
        {"candidate.ndethist": {"$gte": 1}},              # at least one detection in history
        {"properties.sso.is_sso": {"$eq": False}},        # not a solar-system object
        {"candidate.pixelFlags": {"$eq": False}},         # no bad/saturated pixels
        {"candidate.centroid_flag": {"$eq": False}},      # centroid measured cleanly
        {"candidate.extendedness": {"$lte": 0.8}},        # point-like source
        {"candidate.shape_flag": {"$eq": False}},         # shape measured cleanly
        {"candidate.reliability": {"$gt": 0.5}},          # real/bogus reliability score
    ]}},

    # --- 2. Keep only the fields used below ---
    {"$project": {
        "objectId": 1,
        "candidate.jd": 1,
        "candidate.isDipole": 1,
        "candidate.glint_trail": 1,
        "candidate.isdiffpos": 1,
        "properties.rock": 1,
        "properties.stationary": 1,
        "properties.star": 1,
        "properties.near_brightstar": 1,
        "candidate.ndethist": 1,
        "properties.sso.is_sso": 1,
        "candidate.pixelFlags": 1,
        "candidate.centroid_flag": 1,
        "candidate.extendedness": 1,
        "candidate.shape_flag": 1,
        "candidate.reliability": 1,
        "prv_candidates": 1,
        "cross_matches.CatWISE2020": 1,
        "cross_matches.milliquas_v8": 1,
        "cross_matches.Gaia_DR3": 1,
        "cross_matches.LSPSC": 1,
        "cross_matches.PS1_DR2": 1,
    }},

    # --- 3. Cross-match flags ---
    {"$addFields": {
        # Gaia star brighter than G = 16 within 20"
        "any_gaia_brightstar": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.Gaia_DR3", []]},
            "in": {"$and": [
                {"$lte": ["$$this.distance_arcsec", 20]},
                {"$lte": ["$$this.phot_g_mean_mag", 16]},
            ]},
        }}},
        # Star-like LSPSC source (score >= 0.7) close enough for its brightness
        "star_lspsc": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.LSPSC", []]},
            "in": {"$and": [
                {"$gte": ["$$this.score", 0.7]},
                {"$or": [
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 4]},  {"$lte": ["$$this.mag_white", 18]}]},
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 10]}, {"$lte": ["$$this.mag_white", 16]}]},
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 30]}, {"$lte": ["$$this.mag_white", 15]}]},
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 6]},  {"$lte": ["$$this.mag_white", 17.5]}]},
                ]},
            ]},
        }}},
        # Red WISE colour (W1 - W2 > 0.5), typical of AGN
        "WiseW1W2_gt1": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.CatWISE2020", []]},
            "in": {"$gt": [{"$subtract": ["$$this.w1mpro", "$$this.w2mpro"]}, 0.5]},
        }}},
        # Gaia source within 1" with a significant parallax (>= 2 sigma), i.e. a star
        "gaia_star": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.Gaia_DR3", []]},
            "in": {"$and": [
                {"$lte": ["$$this.distance_arcsec", 1]},
                {"$gte": [{"$abs": {"$divide": ["$$this.parallax", "$$this.parallax_error"]}}, 2]},
            ]},
        }}},
        # Time of the earliest previous detection
        "jd_min_prv": {"$min": "$prv_candidates.jd"},
        # PS1 source within 0.5" that is not classified as a galaxy
        "ps1_star": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.PS1_DR2", []]},
            "in": {"$and": [
                {"$lt": ["$$this.distance_arcsec", 0.5]},
                {"$ne": ["$$this.strm_class", "GALAXY"]},
            ]},
        }}},
    }},

    # --- 4. Young, not a star, not an AGN/QSO ---
    {"$match": {"$and": [
        # First detected less than 30 days ago
        {"$expr": {"$lt": [{"$subtract": ["$candidate.jd", "$jd_min_prv"]}, 30]}},
        # No bright Gaia star or star-like LSPSC source nearby
        {"$and": [
            {"any_gaia_brightstar": {"$eq": False}},
            {"star_lspsc": {"$eq": False}},
        ]},
        # No WISE match, or the WISE match is not AGN-coloured
        {"$or": [
            {"$expr": {"$lt": [{"$size": {"$ifNull": ["$cross_matches.CatWISE2020", []]}}, 1]}},
            {"WiseW1W2_gt1": {"$eq": False}},
        ]},
        {"gaia_star": {"$eq": False}},
        {"ps1_star": {"$eq": False}},
        # No Milliquas (QSO/AGN catalog) match
        {"$expr": {"$lt": [{"$size": {"$ifNull": ["$cross_matches.milliquas_v8", []]}}, 1]}},
    ]}},

    # --- 5. Output fields ---
    {"$project": {
        "objectId": 1,
        "candidate.jd": 1,
        "candidate.isDipole": 1,
        "candidate.glint_trail": 1,
        "candidate.isdiffpos": 1,
        "properties.rock": 1,
        "properties.stationary": 1,
        "properties.star": 1,
        "properties.near_brightstar": 1,
        "candidate.ndethist": 1,
        "properties.sso.is_sso": 1,
        "candidate.pixelFlags": 1,
        "candidate.centroid_flag": 1,
        "candidate.extendedness": 1,
        "candidate.shape_flag": 1,
        "candidate.reliability": 1,
        "prv_candidates": 1,
        "cross_matches.CatWISE2020": 1,
        "cross_matches.milliquas_v8": 1,
        "cross_matches.Gaia_DR3": 1,
        "cross_matches.LSPSC": 1,
        "cross_matches.PS1_DR2": 1,
        "any_gaia_brightstar": 1,
        "star_lspsc": 1,
        "WiseW1W2_gt1": 1,
        "gaia_star": 1,
        "ps1_star": 1,
        "jd_min_prv": 1,
    }},
]


# =============================================================================
# LSST_EMGW_BASIC
# Same as LSST_EMGW_NON_NUCLEAR, plus a nuclear veto: transients within 0.5"
# of a galaxy nucleus (NED source, or galaxy-like LSPSC source) are rejected,
# since these are likely AGN/TDEs rather than GW counterparts.
# =============================================================================

LSST_EMGW_BASIC = [
    # --- 1. Alert quality and contaminant cuts ---
    {"$match": {"$and": [
        {"candidate.isDipole": {"$eq": False}},           # not a dipole (bad subtraction)
        {"candidate.glint_trail": {"$eq": False}},        # not a satellite glint / trail
        {"candidate.isdiffpos": {"$eq": True}},           # positive difference flux (brightening)
        {"properties.rock": {"$eq": False}},              # not a known asteroid
        {"properties.stationary": {"$eq": True}},         # not moving between detections
        {"properties.star": {"$eq": False}},              # not a known star
        {"properties.near_brightstar": {"$eq": False}},   # not in a bright-star halo
        {"candidate.ndethist": {"$gte": 1}},              # at least one detection in history
        {"properties.sso.is_sso": {"$eq": False}},        # not a solar-system object
        {"candidate.pixelFlags": {"$eq": False}},         # no bad/saturated pixels
        {"candidate.centroid_flag": {"$eq": False}},      # centroid measured cleanly
        {"candidate.extendedness": {"$lte": 0.8}},        # point-like source
        {"candidate.shape_flag": {"$eq": False}},         # shape measured cleanly
        {"candidate.reliability": {"$gt": 0.5}},          # real/bogus reliability score
    ]}},

    # --- 2. Keep only the fields used below ---
    {"$project": {
        "objectId": 1,
        "candidate.jd": 1,
        "candidate.isDipole": 1,
        "candidate.glint_trail": 1,
        "candidate.isdiffpos": 1,
        "properties.rock": 1,
        "properties.stationary": 1,
        "properties.star": 1,
        "properties.near_brightstar": 1,
        "candidate.ndethist": 1,
        "properties.sso.is_sso": 1,
        "candidate.pixelFlags": 1,
        "candidate.centroid_flag": 1,
        "candidate.extendedness": 1,
        "candidate.shape_flag": 1,
        "candidate.reliability": 1,
        "prv_candidates": 1,
        "cross_matches.CatWISE2020": 1,
        "cross_matches.milliquas_v8": 1,
        "cross_matches.Gaia_DR3": 1,
        "cross_matches.LSPSC": 1,
        "cross_matches.PS1_DR2": 1,
        "cross_matches.NED": 1,
    }},

    # --- 3. Cross-match flags ---
    {"$addFields": {
        # Gaia star brighter than G = 16 within 20"
        "any_gaia_brightstar": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.Gaia_DR3", []]},
            "in": {"$and": [
                {"$lte": ["$$this.distance_arcsec", 20]},
                {"$lte": ["$$this.phot_g_mean_mag", 16]},
            ]},
        }}},
        # Star-like LSPSC source (score >= 0.7) close enough for its brightness
        "star_lspsc": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.LSPSC", []]},
            "in": {"$and": [
                {"$gte": ["$$this.score", 0.7]},
                {"$or": [
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 4]},  {"$lte": ["$$this.mag_white", 18]}]},
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 10]}, {"$lte": ["$$this.mag_white", 16]}]},
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 30]}, {"$lte": ["$$this.mag_white", 15]}]},
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 6]},  {"$lte": ["$$this.mag_white", 17.5]}]},
                ]},
            ]},
        }}},
        # Red WISE colour (W1 - W2 > 0.5), typical of AGN
        "WiseW1W2_gt1": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.CatWISE2020", []]},
            "in": {"$gt": [{"$subtract": ["$$this.w1mpro", "$$this.w2mpro"]}, 0.5]},
        }}},
        # Gaia source within 1" with a significant parallax (>= 2 sigma), i.e. a star
        "gaia_star": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.Gaia_DR3", []]},
            "in": {"$and": [
                {"$lte": ["$$this.distance_arcsec", 1]},
                {"$gte": [{"$abs": {"$divide": ["$$this.parallax", "$$this.parallax_error"]}}, 2]},
            ]},
        }}},
        # Time of the earliest previous detection
        "jd_min_prv": {"$min": "$prv_candidates.jd"},
        # Galaxy-like LSPSC source (score <= 0.5) within 0.5": transient on a galaxy nucleus
        "LSPSCnuclearmatch_copy": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.LSPSC", []]},
            "in": {"$and": [
                {"$lte": ["$$this.distance_arcsec", 0.5]},
                {"$lte": ["$$this.score", 0.5]},
            ]},
        }}},
        # Any NED object within 0.5"
        "NEDnuclearcrossmatch": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.NED", []]},
            "in": {"$lte": ["$$this.distance_arcsec", 0.5]},
        }}},
        # PS1 source within 0.5" that is not classified as a galaxy
        "ps1_star": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.PS1_DR2", []]},
            "in": {"$and": [
                {"$lt": ["$$this.distance_arcsec", 0.5]},
                {"$ne": ["$$this.strm_class", "GALAXY"]},
            ]},
        }}},
    }},

    # --- 4. Young, not a star, not an AGN/QSO, not nuclear ---
    {"$match": {"$and": [
        {"$and": [
            # First detected less than 30 days ago
            {"$expr": {"$lt": [{"$subtract": ["$candidate.jd", "$jd_min_prv"]}, 30]}},
            # No bright Gaia star or star-like LSPSC source nearby
            {"$and": [
                {"any_gaia_brightstar": {"$eq": False}},
                {"star_lspsc": {"$eq": False}},
            ]},
            # No WISE match, or the WISE match is not AGN-coloured
            {"$or": [
                {"$expr": {"$lt": [{"$size": {"$ifNull": ["$cross_matches.CatWISE2020", []]}}, 1]}},
                {"WiseW1W2_gt1": {"$eq": False}},
            ]},
            {"gaia_star": {"$eq": False}},
            {"ps1_star": {"$eq": False}},
            # No Milliquas (QSO/AGN catalog) match
            {"$expr": {"$lt": [{"$size": {"$ifNull": ["$cross_matches.milliquas_v8", []]}}, 1]}},
        ]},
        # Nuclear veto: no NED or galaxy-like LSPSC source within 0.5"
        {"$nor": [{"$or": [
            {"NEDnuclearcrossmatch": {"$eq": True}},
            {"LSPSCnuclearmatch_copy": {"$eq": True}},
        ]}]},
    ]}},

    # --- 5. Output fields ---
    {"$project": {
        "objectId": 1,
        "candidate.jd": 1,
        "candidate.isDipole": 1,
        "candidate.glint_trail": 1,
        "candidate.isdiffpos": 1,
        "properties.rock": 1,
        "properties.stationary": 1,
        "properties.star": 1,
        "properties.near_brightstar": 1,
        "candidate.ndethist": 1,
        "properties.sso.is_sso": 1,
        "candidate.pixelFlags": 1,
        "candidate.centroid_flag": 1,
        "candidate.extendedness": 1,
        "candidate.shape_flag": 1,
        "candidate.reliability": 1,
        "prv_candidates": 1,
        "cross_matches.CatWISE2020": 1,
        "cross_matches.milliquas_v8": 1,
        "cross_matches.Gaia_DR3": 1,
        "cross_matches.LSPSC": 1,
        "cross_matches.PS1_DR2": 1,
        "cross_matches.NED": 1,
        "any_gaia_brightstar": 1,
        "star_lspsc": 1,
        "WiseW1W2_gt1": 1,
        "gaia_star": 1,
        "ps1_star": 1,
        "NEDnuclearcrossmatch": 1,
        "LSPSCnuclearmatch_copy": 1,
        "jd_min_prv": 1,
    }},
]


# =============================================================================
# LSST_EMGW_FAST
# Same as LSST_EMGW_BASIC, and the transient must also be evolving fast, as
# expected for a kilonova. In at least one band (g, r, i, z, y, u), either:
#   - fading  at >= 0.3 mag/day, or
#   - rising  at >= 1   mag/day,
# with the rate significant at >= 2 sigma, measured over >= 0.05 days.
# (Magnitudes: a positive rate is fading, a negative rate is rising.)
#
# NOTE: the u-band fading cut uses "u.fading.dt", while the other bands use
# the band-level "<band>.dt". Kept as it was; check which is intended.
# =============================================================================

LSST_EMGW_FAST = [
    # --- 1. Alert quality and contaminant cuts ---
    {"$match": {"$and": [
        {"candidate.isDipole": {"$eq": False}},           # not a dipole (bad subtraction)
        {"candidate.glint_trail": {"$eq": False}},        # not a satellite glint / trail
        {"candidate.isdiffpos": {"$eq": True}},           # positive difference flux (brightening)
        {"properties.rock": {"$eq": False}},              # not a known asteroid
        {"properties.stationary": {"$eq": True}},         # not moving between detections
        {"properties.star": {"$eq": False}},              # not a known star
        {"properties.near_brightstar": {"$eq": False}},   # not in a bright-star halo
        {"candidate.ndethist": {"$gte": 1}},              # at least one detection in history
        {"properties.sso.is_sso": {"$eq": False}},        # not a solar-system object
        {"candidate.pixelFlags": {"$eq": False}},         # no bad/saturated pixels
        {"candidate.centroid_flag": {"$eq": False}},      # centroid measured cleanly
        {"candidate.extendedness": {"$lte": 0.8}},        # point-like source
        {"candidate.shape_flag": {"$eq": False}},         # shape measured cleanly
        {"candidate.reliability": {"$gt": 0.5}},          # real/bogus reliability score
    ]}},

    # --- 2. Keep only the fields used below ---
    {"$project": {
        "objectId": 1,
        "candidate.jd": 1,
        "candidate.isDipole": 1,
        "candidate.glint_trail": 1,
        "candidate.isdiffpos": 1,
        "properties.rock": 1,
        "properties.stationary": 1,
        "properties.star": 1,
        "properties.near_brightstar": 1,
        "candidate.ndethist": 1,
        "properties.sso.is_sso": 1,
        "candidate.pixelFlags": 1,
        "candidate.centroid_flag": 1,
        "candidate.extendedness": 1,
        "candidate.shape_flag": 1,
        "candidate.reliability": 1,
        "prv_candidates": 1,
        "cross_matches.CatWISE2020": 1,
        "cross_matches.milliquas_v8": 1,
        "cross_matches.Gaia_DR3": 1,
        "cross_matches.LSPSC": 1,
        "cross_matches.PS1_DR2": 1,
        "cross_matches.NED": 1,
        # light-curve rates per band
        "properties.multisurvey_photstats.g.fading.rate": 1,
        "properties.multisurvey_photstats.g.fading.rate_error": 1,
        "properties.multisurvey_photstats.g.dt": 1,
        "properties.multisurvey_photstats.g.rising.rate": 1,
        "properties.multisurvey_photstats.g.rising.rate_error": 1,
        "properties.multisurvey_photstats.g.rising.dt": 1,
        "properties.multisurvey_photstats.r.fading.rate": 1,
        "properties.multisurvey_photstats.r.fading.rate_error": 1,
        "properties.multisurvey_photstats.r.dt": 1,
        "properties.multisurvey_photstats.r.rising.rate": 1,
        "properties.multisurvey_photstats.r.rising.rate_error": 1,
        "properties.multisurvey_photstats.r.rising.dt": 1,
        "properties.multisurvey_photstats.i.fading.rate": 1,
        "properties.multisurvey_photstats.i.fading.rate_error": 1,
        "properties.multisurvey_photstats.i.dt": 1,
        "properties.multisurvey_photstats.i.rising.rate": 1,
        "properties.multisurvey_photstats.i.rising.rate_error": 1,
        "properties.multisurvey_photstats.i.rising.dt": 1,
        "properties.multisurvey_photstats.z.fading.rate": 1,
        "properties.multisurvey_photstats.z.fading.rate_error": 1,
        "properties.multisurvey_photstats.z.dt": 1,
        "properties.multisurvey_photstats.z.rising.rate": 1,
        "properties.multisurvey_photstats.z.rising.rate_error": 1,
        "properties.multisurvey_photstats.z.rising.dt": 1,
        "properties.multisurvey_photstats.y.fading.rate": 1,
        "properties.multisurvey_photstats.y.fading.rate_error": 1,
        "properties.multisurvey_photstats.y.dt": 1,
        "properties.multisurvey_photstats.y.rising.rate": 1,
        "properties.multisurvey_photstats.y.rising.rate_error": 1,
        "properties.multisurvey_photstats.y.rising.dt": 1,
        "properties.multisurvey_photstats.u.fading.rate": 1,
        "properties.multisurvey_photstats.u.fading.rate_error": 1,
        "properties.multisurvey_photstats.u.fading.dt": 1,
        "properties.multisurvey_photstats.u.rising.rate": 1,
        "properties.multisurvey_photstats.u.rising.rate_error": 1,
        "properties.multisurvey_photstats.u.rising.dt": 1,
    }},

    # --- 3. Cross-match flags ---
    {"$addFields": {
        # Gaia star brighter than G = 16 within 20"
        "any_gaia_brightstar": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.Gaia_DR3", []]},
            "in": {"$and": [
                {"$lte": ["$$this.distance_arcsec", 20]},
                {"$lte": ["$$this.phot_g_mean_mag", 16]},
            ]},
        }}},
        # Star-like LSPSC source (score >= 0.7) close enough for its brightness
        "star_lspsc": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.LSPSC", []]},
            "in": {"$and": [
                {"$gte": ["$$this.score", 0.7]},
                {"$or": [
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 4]},  {"$lte": ["$$this.mag_white", 18]}]},
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 10]}, {"$lte": ["$$this.mag_white", 16]}]},
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 30]}, {"$lte": ["$$this.mag_white", 15]}]},
                    {"$and": [{"$lte": ["$$this.distance_arcsec", 6]},  {"$lte": ["$$this.mag_white", 17.5]}]},
                ]},
            ]},
        }}},
        # Red WISE colour (W1 - W2 > 0.5), typical of AGN
        "WiseW1W2_gt1": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.CatWISE2020", []]},
            "in": {"$gt": [{"$subtract": ["$$this.w1mpro", "$$this.w2mpro"]}, 0.5]},
        }}},
        # Gaia source within 1" with a significant parallax (>= 2 sigma), i.e. a star
        "gaia_star": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.Gaia_DR3", []]},
            "in": {"$and": [
                {"$lte": ["$$this.distance_arcsec", 1]},
                {"$gte": [{"$abs": {"$divide": ["$$this.parallax", "$$this.parallax_error"]}}, 2]},
            ]},
        }}},
        # Time of the earliest previous detection
        "jd_min_prv": {"$min": "$prv_candidates.jd"},
        # Galaxy-like LSPSC source (score <= 0.5) within 0.5": transient on a galaxy nucleus
        "LSPSCnuclearmatch_copy": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.LSPSC", []]},
            "in": {"$and": [
                {"$lte": ["$$this.distance_arcsec", 0.5]},
                {"$lte": ["$$this.score", 0.5]},
            ]},
        }}},
        # Any NED object within 0.5"
        "NEDnuclearcrossmatch": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.NED", []]},
            "in": {"$lte": ["$$this.distance_arcsec", 0.5]},
        }}},
        # PS1 source within 0.5" that is not classified as a galaxy
        "ps1_star": {"$anyElementTrue": {"$map": {
            "input": {"$ifNull": ["$cross_matches.PS1_DR2", []]},
            "in": {"$and": [
                {"$lt": ["$$this.distance_arcsec", 0.5]},
                {"$ne": ["$$this.strm_class", "GALAXY"]},
            ]},
        }}},
    }},

    # --- 4. Young, not a star, not an AGN/QSO, not nuclear, evolving fast ---
    {"$match": {"$and": [
        {"$and": [
            # First detected less than 30 days ago
            {"$expr": {"$lt": [{"$subtract": ["$candidate.jd", "$jd_min_prv"]}, 30]}},
            # No bright Gaia star or star-like LSPSC source nearby
            {"$and": [
                {"any_gaia_brightstar": {"$eq": False}},
                {"star_lspsc": {"$eq": False}},
            ]},
            # No WISE match, or the WISE match is not AGN-coloured
            {"$or": [
                {"$expr": {"$lt": [{"$size": {"$ifNull": ["$cross_matches.CatWISE2020", []]}}, 1]}},
                {"WiseW1W2_gt1": {"$eq": False}},
            ]},
            {"gaia_star": {"$eq": False}},
            {"ps1_star": {"$eq": False}},
            # No Milliquas (QSO/AGN catalog) match
            {"$expr": {"$lt": [{"$size": {"$ifNull": ["$cross_matches.milliquas_v8", []]}}, 1]}},
        ]},
        # Nuclear veto: no NED or galaxy-like LSPSC source within 0.5"
        {"$nor": [{"$or": [
            {"NEDnuclearcrossmatch": {"$eq": True}},
            {"LSPSCnuclearmatch_copy": {"$eq": True}},
        ]}]},
        # Fast evolution in at least one band
        {"$or": [
            # Fading >= 0.3 mag/day at >= 2 sigma
            {"$and": [
                {"properties.multisurvey_photstats.g.fading.rate": {"$gte": 0.3}},
                {"properties.multisurvey_photstats.g.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    "$properties.multisurvey_photstats.g.fading.rate",
                    "$properties.multisurvey_photstats.g.fading.rate_error",
                ]}, 2]}},
            ]},
            {"$and": [
                {"properties.multisurvey_photstats.r.fading.rate": {"$gte": 0.3}},
                {"properties.multisurvey_photstats.r.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    "$properties.multisurvey_photstats.r.fading.rate",
                    "$properties.multisurvey_photstats.r.fading.rate_error",
                ]}, 2]}},
            ]},
            {"$and": [
                {"properties.multisurvey_photstats.i.fading.rate": {"$gte": 0.3}},
                {"properties.multisurvey_photstats.i.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    "$properties.multisurvey_photstats.i.fading.rate",
                    "$properties.multisurvey_photstats.i.fading.rate_error",
                ]}, 2]}},
            ]},
            {"$and": [
                {"properties.multisurvey_photstats.z.fading.rate": {"$gte": 0.3}},
                {"properties.multisurvey_photstats.z.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    "$properties.multisurvey_photstats.z.fading.rate",
                    "$properties.multisurvey_photstats.z.fading.rate_error",
                ]}, 2]}},
            ]},
            {"$and": [
                {"properties.multisurvey_photstats.y.fading.rate": {"$gte": 0.3}},
                {"properties.multisurvey_photstats.y.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    "$properties.multisurvey_photstats.y.fading.rate",
                    "$properties.multisurvey_photstats.y.fading.rate_error",
                ]}, 2]}},
            ]},
            {"$and": [
                {"properties.multisurvey_photstats.u.fading.rate": {"$gte": 0.3}},
                {"properties.multisurvey_photstats.u.fading.dt": {"$gte": 0.05}},  # see NOTE above
                {"$expr": {"$gte": [{"$divide": [
                    "$properties.multisurvey_photstats.u.fading.rate",
                    "$properties.multisurvey_photstats.u.fading.rate_error",
                ]}, 2]}},
            ]},
            # Rising >= 1 mag/day at >= 2 sigma
            {"$and": [
                {"properties.multisurvey_photstats.g.rising.rate": {"$lte": -1}},
                {"properties.multisurvey_photstats.g.rising.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    {"$abs": "$properties.multisurvey_photstats.g.rising.rate"},
                    "$properties.multisurvey_photstats.g.rising.rate_error",
                ]}, 2]}},
            ]},
            {"$and": [
                {"properties.multisurvey_photstats.r.rising.rate": {"$lte": -1}},
                {"properties.multisurvey_photstats.r.rising.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    {"$abs": "$properties.multisurvey_photstats.r.rising.rate"},
                    "$properties.multisurvey_photstats.r.rising.rate_error",
                ]}, 2]}},
            ]},
            {"$and": [
                {"properties.multisurvey_photstats.i.rising.rate": {"$lte": -1}},
                {"properties.multisurvey_photstats.i.rising.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    {"$abs": "$properties.multisurvey_photstats.i.rising.rate"},
                    "$properties.multisurvey_photstats.i.rising.rate_error",
                ]}, 2]}},
            ]},
            {"$and": [
                {"properties.multisurvey_photstats.z.rising.rate": {"$lte": -1}},
                {"properties.multisurvey_photstats.z.rising.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    {"$abs": "$properties.multisurvey_photstats.z.rising.rate"},
                    "$properties.multisurvey_photstats.z.rising.rate_error",
                ]}, 2]}},
            ]},
            {"$and": [
                {"properties.multisurvey_photstats.y.rising.rate": {"$lte": -1}},
                {"properties.multisurvey_photstats.y.rising.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    {"$abs": "$properties.multisurvey_photstats.y.rising.rate"},
                    "$properties.multisurvey_photstats.y.rising.rate_error",
                ]}, 2]}},
            ]},
            {"$and": [
                {"properties.multisurvey_photstats.u.rising.rate": {"$lte": -1}},
                {"properties.multisurvey_photstats.u.rising.dt": {"$gte": 0.05}},
                {"$expr": {"$gte": [{"$divide": [
                    {"$abs": "$properties.multisurvey_photstats.u.rising.rate"},
                    "$properties.multisurvey_photstats.u.rising.rate_error",
                ]}, 2]}},
            ]},
        ]},
    ]}},

    # --- 5. Output fields ---
    {"$project": {
        "objectId": 1,
        "candidate.jd": 1,
        "candidate.isDipole": 1,
        "candidate.glint_trail": 1,
        "candidate.isdiffpos": 1,
        "properties.rock": 1,
        "properties.stationary": 1,
        "properties.star": 1,
        "properties.near_brightstar": 1,
        "candidate.ndethist": 1,
        "properties.sso.is_sso": 1,
        "candidate.pixelFlags": 1,
        "candidate.centroid_flag": 1,
        "candidate.extendedness": 1,
        "candidate.shape_flag": 1,
        "candidate.reliability": 1,
        "prv_candidates": 1,
        "cross_matches.CatWISE2020": 1,
        "cross_matches.milliquas_v8": 1,
        "cross_matches.Gaia_DR3": 1,
        "cross_matches.LSPSC": 1,
        "cross_matches.PS1_DR2": 1,
        "cross_matches.NED": 1,
        # light-curve rates per band
        "properties.multisurvey_photstats.g.fading.rate": 1,
        "properties.multisurvey_photstats.g.fading.rate_error": 1,
        "properties.multisurvey_photstats.g.dt": 1,
        "properties.multisurvey_photstats.g.rising.rate": 1,
        "properties.multisurvey_photstats.g.rising.rate_error": 1,
        "properties.multisurvey_photstats.g.rising.dt": 1,
        "properties.multisurvey_photstats.r.fading.rate": 1,
        "properties.multisurvey_photstats.r.fading.rate_error": 1,
        "properties.multisurvey_photstats.r.dt": 1,
        "properties.multisurvey_photstats.r.rising.rate": 1,
        "properties.multisurvey_photstats.r.rising.rate_error": 1,
        "properties.multisurvey_photstats.r.rising.dt": 1,
        "properties.multisurvey_photstats.i.fading.rate": 1,
        "properties.multisurvey_photstats.i.fading.rate_error": 1,
        "properties.multisurvey_photstats.i.dt": 1,
        "properties.multisurvey_photstats.i.rising.rate": 1,
        "properties.multisurvey_photstats.i.rising.rate_error": 1,
        "properties.multisurvey_photstats.i.rising.dt": 1,
        "properties.multisurvey_photstats.z.fading.rate": 1,
        "properties.multisurvey_photstats.z.fading.rate_error": 1,
        "properties.multisurvey_photstats.z.dt": 1,
        "properties.multisurvey_photstats.z.rising.rate": 1,
        "properties.multisurvey_photstats.z.rising.rate_error": 1,
        "properties.multisurvey_photstats.z.rising.dt": 1,
        "properties.multisurvey_photstats.y.fading.rate": 1,
        "properties.multisurvey_photstats.y.fading.rate_error": 1,
        "properties.multisurvey_photstats.y.dt": 1,
        "properties.multisurvey_photstats.y.rising.rate": 1,
        "properties.multisurvey_photstats.y.rising.rate_error": 1,
        "properties.multisurvey_photstats.y.rising.dt": 1,
        "properties.multisurvey_photstats.u.fading.rate": 1,
        "properties.multisurvey_photstats.u.fading.rate_error": 1,
        "properties.multisurvey_photstats.u.fading.dt": 1,
        "properties.multisurvey_photstats.u.rising.rate": 1,
        "properties.multisurvey_photstats.u.rising.rate_error": 1,
        "properties.multisurvey_photstats.u.rising.dt": 1,
        # computed flags
        "any_gaia_brightstar": 1,
        "star_lspsc": 1,
        "WiseW1W2_gt1": 1,
        "gaia_star": 1,
        "ps1_star": 1,
        "NEDnuclearcrossmatch": 1,
        "LSPSCnuclearmatch_copy": 1,
        "jd_min_prv": 1,
    }},
]


# =============================================================================
# Registry
# =============================================================================

all_filters = {
    "ZTF_BASIC": ZTF_BASIC,
    "LSST_EMGW_NON_NUCLEAR": LSST_EMGW_NON_NUCLEAR,
    "LSST_EMGW_BASIC": LSST_EMGW_BASIC,
    "LSST_EMGW_FAST": LSST_EMGW_FAST,
}

# fritz filter ids for these filters (helps in updating filters)
filter_ids = {
    "ZTF_BASIC":1422, 
    "LSST_EMGW_NON_NUCLEAR":1441, 
    "LSST_EMGW_BASIC":1442, 
    "LSST_EMGW_FAST":1226
}



def list_filters():
    return list(all_filters.keys())
