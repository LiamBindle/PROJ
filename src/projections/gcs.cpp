#define PJ_LIB__

#include <errno.h>
#include <math.h>

#include "proj.h"
#include "proj_internal.h"
#include <math.h>

PROJ_HEAD(gcs, "GnominicCubedSphere") "\n\tAzi, Sph";

#define EPS10  1.e-10

PJ_XY decode_tile_xy(char c, const char** error_msg);
void pcns1_tile_rotations(PJ_XY* tile_pos, int* ccw_rotations, const char** error_msg);
PJ *gcs_destructor(PJ *P, int errlev);

static inline PJ_XYZ sph2cart(const PJ_LP& lp);
static inline PJ_LP cart2sph(const PJ_XYZ& xyz);
static inline void rotate_xyz(PJ_XYZ& xyz, const double rot[3][3], bool fwd);

inline void fwd_schmidt(PJ_LP& lp, double s);
inline void inv_schmidt(PJ_LP& lp, double s);

inline void fwd_schmidt(PJ_LP& lp, double s);
inline void inv_schmidt(PJ_LP& lp, double s);

namespace { // anonymous namespace

// xyz->uv transformation matrix for the base net (bn); a base net is one that requires no rotations
constexpr double bn_trans_mat_xyz_to_uv[6][2][3] = { 
    {{0, 1, 0},  {0, 0, 1}},  // +X
    {{0, -1, 0}, {0, 0, 1}},  // -X
    {{-1, 0, 0}, {0, 0, 1}},  // +Y
    {{1, 0, 0},  {0, 0, 1}},  // -Y
    {{0, 1, 0},  {-1, 0, 0}}, // +Z
    {{0, 1, 0},  {1, 0, 0}}   // -Z
};

constexpr double face_xyz[6][3] = {
    {1, 0, 0},  // +X
    {-1, 0, 0}, // -X
    {0, 1, 0},  // +Y
    {0, -1, 0}, // -Y
    {0, 0, 1},  // +Z
    {0, 0, -1}  // -Z
};

enum Face {
    X_P = 0,
    X_N = 1,
    Y_P = 2,
    Y_N = 3,
    Z_P = 4,
    Z_N = 5,
    NOT_A_FACE=-1
};
} // anonymous namespace

namespace { // anonymous namespace
struct pj_opaque {
    bool specific_face;                 // user selected a specific face
    enum Face selected_face;            // the specific face selected by the user
    PJ_XY tile_false_origin[6];         // false origin of each tile (XY space)
    int tile_rotation_mat[6][2][2];     // forward rotation matrix for each face
    PJ_LP face_false_origin[6];         // false origin of each face (LP space)
    PJ* face_proj[6];                   // gnom proj object for each face
    char face_proj_def[6][256];         // gnom proj definition string of each face
    enum Face mru_face[6];              // most recently used face

    double stretch_factor;
    double fwd_rot_z[3][3];
    double fwd_rot_y[3][3];
    bool target_pt_rot;
};
} // anonymous namespace

PJ_XY decode_tile_xy(char c, const char** error_msg) {
    PJ_XY xy = {0.0,0.0};
    int v = toupper(c) - 'A';
    *error_msg = nullptr;

    if (v < 0 || v > 18) {
        *error_msg = "Invalid character in +cube_net (valid characters are A-S case-insensitive)";
        return xy;
    }
    
    // XY offset
    if (v < 10) {
        xy.x = static_cast<double>(v % 5);
        xy.y = -static_cast<double>(v / 5);
    } else if (v < 14) {
        xy.x = static_cast<double>(v - 10);
        xy.y = -2;
    } else if (v < 17) {
        xy.x = static_cast<double>(v - 14);
        xy.y = -3;
    } else {
        xy.x = static_cast<double>(v - 17);
        xy.y = -4;
    }

    // multiply by 2
    xy.x *= 2.;
    xy.y *= 2.;
    return xy;
}

void pcns1_tile_rotations(PJ_XY* tile_pos, int* ccw_rotations, const char** error_msg) {
    *error_msg = nullptr;
    // calculate rotations for each face
    for (int i = 0; i < 6; ++i) {
        int j;
        PJ_XYZ anchor={0.,0.,0.};
        PJ_UV base_edge={0.,0.}, target_edge={0.,0.};
        double l1norm;
        double temp;

        // find an adjacent tile (it doesn't matter which)
        for (j = 0; j < 6; ++j) {
            l1norm = fabs(tile_pos[i].x - tile_pos[j].x) + fabs(tile_pos[i].y - tile_pos[j].y);
            if(fabs(l1norm - 2.) < EPS10) {
                target_edge.u = (tile_pos[i].x + tile_pos[j].x) / 2. - tile_pos[i].x;
                target_edge.v = (tile_pos[i].y + tile_pos[j].y) / 2. - tile_pos[i].y;
                break; /* tiles i,j are adjacent */
            }
        }
        if (j == 6) {
            *error_msg = "Invalid +cube_net (one or more tiles has no adjacent tiles)";
            return;
        }

        // calculate the coordinate where cube face i and j connect
        anchor.x = face_xyz[i][0] + face_xyz[j][0];
        anchor.y = face_xyz[i][1] + face_xyz[j][1];
        anchor.z = face_xyz[i][2] + face_xyz[j][2];

        // determine edge of face i in tile-centered coords (UV) that connects face j
        base_edge.u = bn_trans_mat_xyz_to_uv[i][0][0] * (anchor.x - face_xyz[i][0])
                      + bn_trans_mat_xyz_to_uv[i][0][1] * (anchor.y - face_xyz[i][1])
                      + bn_trans_mat_xyz_to_uv[i][0][2] * (anchor.z - face_xyz[i][2]);
        base_edge.v = bn_trans_mat_xyz_to_uv[i][1][0] * (anchor.x - face_xyz[i][0])
                      + bn_trans_mat_xyz_to_uv[i][1][1] * (anchor.y - face_xyz[i][1])
                      + bn_trans_mat_xyz_to_uv[i][1][2] * (anchor.z - face_xyz[i][2]);

        ccw_rotations[i] = 0;
        l1norm = fabs(base_edge.u - target_edge.u) + fabs(base_edge.v - target_edge.v);
        while(l1norm > EPS10) { /* native edge is not equal to target edge */
            /* rotate 90 deg CCW */
            temp = -base_edge.v;
            base_edge.v = base_edge.u;
            base_edge.u = temp;
            ccw_rotations[i] += 1;

            l1norm = fabs(base_edge.u - target_edge.u) + fabs(base_edge.v - target_edge.v);

            if(ccw_rotations[i] > 3) {
                *error_msg = "+cube_net is not a valid cube net";
                return;
            }
        }
    }

}

static inline PJ_XYZ sph2cart(const PJ_LP& lp) {
    PJ_XYZ xyz;
    xyz.x = cos(lp.phi) * cos(lp.lam);
    xyz.y = cos(lp.phi) * sin(lp.lam);
    xyz.z =  sin(lp.phi);
    return xyz;
}

static inline PJ_LP cart2sph(const PJ_XYZ& xyz) {
    PJ_LP lp;
    lp.phi = asin(xyz.z);
    lp.lam = atan2(xyz.y, xyz.x);
    return lp;
}

static inline void rotate_xyz(PJ_XYZ& xyz, const double rot[3][3], bool fwd=true) {
    PJ_XYZ temp_xyz;

    if(fwd) {
        temp_xyz.x = rot[0][0]*xyz.x + rot[0][1]*xyz.y + rot[0][2]*xyz.z;
        temp_xyz.y = rot[1][0]*xyz.x + rot[1][1]*xyz.y + rot[1][2]*xyz.z;
        temp_xyz.z = rot[2][0]*xyz.x + rot[2][1]*xyz.y + rot[2][2]*xyz.z;
    } else {
        temp_xyz.x = rot[0][0]*xyz.x + rot[1][0]*xyz.y + rot[2][0]*xyz.z;
        temp_xyz.y = rot[0][1]*xyz.x + rot[1][1]*xyz.y + rot[2][1]*xyz.z;
        temp_xyz.z = rot[0][2]*xyz.x + rot[1][2]*xyz.y + rot[2][2]*xyz.z;
    }

    xyz.x = temp_xyz.x;
    xyz.y = temp_xyz.y;
    xyz.z = temp_xyz.z;
}

inline void fwd_schmidt(PJ_LP& lp, double s) {
    double D = (1. - s*s) / (1. + s*s);
    lp.phi = asin((D + sin(lp.phi)) / (1 + D * sin(lp.phi)));
}

inline void inv_schmidt(PJ_LP& lp, double s) {
    double D = (1. - s*s) / (1. + s*s);
    lp.phi = -asin((sin(lp.phi) - D) / (D * sin(lp.phi) - 1.));
}

static PJ_XY gcs_s_forward (PJ_LP lp, PJ *P) {           /* Spheroidal, forward */
    PJ_XY xy = {0.0,0.0};
    struct pj_opaque *Q = static_cast<struct pj_opaque*>(P->opaque);
    double temp;
    enum Face i;
    PJ_LP lp_i;
    int proj_errno;

    if (Q->target_pt_rot) {
        PJ_XYZ xyz = sph2cart(lp);
        rotate_xyz(xyz, Q->fwd_rot_z);
        rotate_xyz(xyz, Q->fwd_rot_y);
        lp = cart2sph(xyz);
    }

    if (Q->stretch_factor > 1.0) {
        inv_schmidt(lp, Q->stretch_factor);
    }

    int s_end = Q->selected_face ? 1 : 6;  // if a face is selected it is first in MRU list, t.f., stop search for s>=1

    for(int s = 0; s < s_end; ++s) {
        i = static_cast<enum Face>(Q->mru_face[s]);
        lp_i = {lp.lam - Q->face_false_origin[i].lam, lp.phi};
        xy = Q->face_proj[i]->fwd(lp_i, Q->face_proj[i]);
        // correct face if proj errno is 0, and -1 <= x <= 1, and -1 <= y <= 1
        proj_errno = proj_errno_reset(Q->face_proj[i]);
        if (proj_errno == 0 && xy.x >= -1 && xy.x <= 1 && xy.y >= -1 && xy.y <= 1) {
            // Update MRU list
            for(; s >= 1; --s) {
                Q->mru_face[s] = Q->mru_face[s-1];
            }
            Q->mru_face[0] = i;
            break;
        }
    }

    if (proj_errno != 0 || xy.x > 1 || xy.x < -1 || xy.y > 1 || xy.y < -1) {
        proj_errno_set(P, PROJ_ERR_COORD_TRANSFM_OUTSIDE_PROJECTION_DOMAIN);
        xy = {HUGE_VAL, HUGE_VAL};
        return xy;
    }

    /* apply CCW rotations */
    temp = Q->tile_rotation_mat[i][0][0] * xy.x + Q->tile_rotation_mat[i][0][1] * xy.y;
    xy.y = Q->tile_rotation_mat[i][1][0] * xy.x + Q->tile_rotation_mat[i][1][1] * xy.y;
    xy.x = temp;

    /* apply offset */
    xy.x += Q->tile_false_origin[i].x;
    xy.y += Q->tile_false_origin[i].y;

    return xy;
}


static PJ_LP gcs_s_inverse (PJ_XY xy, PJ *P) {           /* Spheroidal, inverse */
    PJ_LP lp = {0.0,0.0};
    struct pj_opaque *Q = static_cast<struct pj_opaque*>(P->opaque);
    double temp;
    double temp_x, temp_y;
    int i;

    for(i = X_P; i <= Z_N; ++i) {
        /* reverse offset */
        temp_x = xy.x - Q->tile_false_origin[i].x;
        temp_y = xy.y - Q->tile_false_origin[i].y;
        if (temp_x >= -1 && temp_x <= 1 && temp_y >= -1 && temp_y <= 1) {
            break;
        }
    }

    if (temp_x > 1 || temp_x < -1 || temp_y > 1 || temp_y < -1) {
        proj_errno_set(P, PROJ_ERR_COORD_TRANSFM_OUTSIDE_PROJECTION_DOMAIN);
        lp = {HUGE_VAL, HUGE_VAL};
        return lp;
    }

    /* reverse rotation */
    temp   = Q->tile_rotation_mat[i][0][0] * temp_x + Q->tile_rotation_mat[i][1][0] * temp_y;
    temp_y = Q->tile_rotation_mat[i][0][1] * temp_x + Q->tile_rotation_mat[i][1][1] * temp_y;
    temp_x = temp;

    lp = Q->face_proj[i]->inv({temp_x, temp_y}, Q->face_proj[i]);
    lp.lam += Q->face_false_origin[i].lam;

    if (Q->stretch_factor > 1.0) {
        fwd_schmidt(lp, Q->stretch_factor);
    }

    if (Q->target_pt_rot) {
        PJ_XYZ xyz = sph2cart(lp);
        rotate_xyz(xyz, Q->fwd_rot_y, false);
        rotate_xyz(xyz, Q->fwd_rot_z, false);
        lp = cart2sph(xyz);
    }
    return lp;
}

PJ *gcs_destructor(PJ *P, int errlev) {
    if (P == nullptr)
        return nullptr;

    if (P->opaque == nullptr)
        return pj_default_destructor (P, errlev);

    struct pj_opaque *Q = static_cast<struct pj_opaque*>(P->opaque);
    for(int i = X_P; i <= Z_N; ++i) {
        if (Q->face_proj[i] != nullptr) {
            Q->face_proj[i]->destructor(Q->face_proj[i], errlev);
        }
    }

    return pj_default_destructor(P, errlev);
}

PJ *PROJECTION(gcs) {
    struct pj_opaque *Q = static_cast<struct pj_opaque*>(calloc (1, sizeof (struct pj_opaque)));
    enum Face cube_index[7] = {NOT_A_FACE, NOT_A_FACE, NOT_A_FACE, NOT_A_FACE, NOT_A_FACE, NOT_A_FACE, NOT_A_FACE}; // mapping user index -> actual index
    char cube_net[6] = {'K', 'H', 'L', 'C', 'G', 'D'};                    // cube net definition
    const char* err_msg = nullptr;
    if (nullptr==Q) {
        return pj_default_destructor (P, PROJ_ERR_OTHER /*ENOMEM*/);
    }
    P->opaque = Q;

    for (int i = 0; i < 6; ++i) {
        Q->face_proj[i] = nullptr; // initialize to null so destructor knows whether to deallocate
    }

    // Get the user's face indices
    if (pj_param(P->ctx, P->params, "tcube_index").i) {
        const char* arg_cube_index = pj_param(P->ctx, P->params, "scube_index").s;
        int not_a_face_count, temp;

        if (strlen(arg_cube_index) != 6) {
            proj_log_error(P, _("+cube_index must be exactly 6 characters"));
            return gcs_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }

        for(int i = X_P; i <= Z_N; ++i) {
            temp = static_cast<int>(arg_cube_index[i] - '0');
            if (temp < 0 || temp > 6) {
                proj_log_error(P, _("Invalid index in +cube_index. Indices must be 0-6."));
                return gcs_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
            }
            cube_index[temp] = static_cast<enum Face>(i);
        }

        // Assert that indices are unique (count of cube_index NOT_A_FACE must be 1)
        not_a_face_count = 0;
        for (int i = 0; i < 7; ++i) {
            if (cube_index[i] == NOT_A_FACE) {
                not_a_face_count += 1;
            }
        }
        if (not_a_face_count != 1) {
            proj_log_error(P, _("Invalid indices in +cube_index. Indices must be unique."));
            return gcs_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }
    } else {
        cube_index[0] = X_P;
        cube_index[1] = Y_P;
        cube_index[2] = Z_P;
        cube_index[3] = X_N;
        cube_index[4] = Y_N;
        cube_index[5] = Z_N;
    }

    // Set specific face if the user specified +cube_face
    Q->specific_face = pj_param(P->ctx, P->params, "tcube_face").i;
    if (Q->specific_face) {
        int temp = pj_param(P->ctx, P->params, "icube_face").i;
        if(temp < 0 || temp >= 7 || cube_index[temp] == NOT_A_FACE) {
            proj_log_error(P, _("Invalid +cube_face (index doesn't occur in +cube_index)."));
            return gcs_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }
        Q->selected_face =  cube_index[temp];
    }

    // Get the cube net definition
    if (pj_param(P->ctx, P->params, "tcube_net").i) {
        const char* temp_cc = pj_param(P->ctx, P->params, "scube_net").s;

        if (strlen(temp_cc) != 6) {
            proj_log_error(P, _("+cube_net must be exactly 6 characters"));
            return gcs_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }

        for (int i = X_P; i <= Z_N; ++i) {
            cube_net[i] = temp_cc[i];
        }
    }
    // Decode the tile positions (in XY space)
    for (int i = X_P; i <= Z_N; ++i) {
        Q->tile_false_origin[i] = decode_tile_xy(cube_net[i], &err_msg);
        if (err_msg != nullptr) {
            proj_log_error(P, _(err_msg));
            return gcs_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }
    }
    // Set the tiles' false origins to the center of +X
    for (int i = Z_N; i >= X_P; --i) {
        Q->tile_false_origin[i].x -= Q->tile_false_origin[X_P].x;
        Q->tile_false_origin[i].y -= Q->tile_false_origin[X_P].y;
    }
    // Solve the tile rotations (number of CCW rotations for each tile)
    int ccw_rotations[6];
    pcns1_tile_rotations(Q->tile_false_origin, ccw_rotations, &err_msg);
    if (err_msg != nullptr) {
        proj_log_error(P, _(err_msg));
        return gcs_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
    }
    // Set tile rotations matrices
    for (int i = X_P; i <= Z_N; ++i) {
        switch (ccw_rotations[i]) {
            case 0:
                Q->tile_rotation_mat[i][0][0] = 1;
                Q->tile_rotation_mat[i][0][1] = 0;
                Q->tile_rotation_mat[i][1][0] = 0;
                Q->tile_rotation_mat[i][1][1] = 1;
                break;
            case 1:
                Q->tile_rotation_mat[i][0][0] = 0;
                Q->tile_rotation_mat[i][0][1] = -1;
                Q->tile_rotation_mat[i][1][0] = 1;
                Q->tile_rotation_mat[i][1][1] = 0;
                break;
            case 2:
                Q->tile_rotation_mat[i][0][0] = -1;
                Q->tile_rotation_mat[i][0][1] = 0;
                Q->tile_rotation_mat[i][1][0] = 0;
                Q->tile_rotation_mat[i][1][1] = -1;
                break;
            case 3:
                Q->tile_rotation_mat[i][0][0] = 0;
                Q->tile_rotation_mat[i][0][1] = 1;
                Q->tile_rotation_mat[i][1][0] = -1;
                Q->tile_rotation_mat[i][1][1] = 0;
                break;
        }
    }

    // Create gnominic projection object for each face
    Q->face_false_origin[X_P] = {0., 0.};
    Q->face_false_origin[X_N] = {M_PI, 0.};
    Q->face_false_origin[Y_P] = {M_HALFPI, 0.};
    Q->face_false_origin[Y_N] = {-M_HALFPI, 0.};
    Q->face_false_origin[Z_P] = {0., M_HALFPI};
    Q->face_false_origin[Z_N] = {0., -M_HALFPI};
    for (int i = X_P; i <= Z_N; ++i) {
        double phi0 = Q->face_false_origin[i].phi*RAD_TO_DEG;
        double lam0 = Q->face_false_origin[i].lam*RAD_TO_DEG;
        sprintf(Q->face_proj_def[i], "+proj=gnom +lat_0=%f +lon_0=%f", phi0, lam0);
        Q->face_proj[i] = proj_create(P->ctx, Q->face_proj_def[i]);
        if (Q->face_proj[i] == nullptr) {
            return gcs_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }
    }

        // Set the tiles' false origins to the center of +X          +cube_face=2 +center_origin
    if (pj_param(P->ctx, P->params, "tcenter_origin").i) {
        if (!pj_param(P->ctx, P->params, "tcube_face").i) {
            proj_log_error(P, _("+center_origin is only valid if +cube_face was specified (no face to center on)"));
            return gcs_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }
        for (int i = X_P; i <= Z_N; ++i) {
            if (i == Q->selected_face) continue;
            Q->tile_false_origin[i].x -= Q->tile_false_origin[Q->selected_face].x;
            Q->tile_false_origin[i].y -= Q->tile_false_origin[Q->selected_face].y;
        }
        Q->tile_false_origin[Q->selected_face].x = 0;
        Q->tile_false_origin[Q->selected_face].y = 0;
    }

    // Initialize MRU list
    for (int i = X_P; i <= Z_N; ++i) {
        Q->mru_face[i] = Q->specific_face ? Q->selected_face : static_cast<enum Face>(i);
    }

    Q->target_pt_rot = false;
    if (pj_param(P->ctx, P->params, "ttlon").i) {
        double theta_z = -pj_param(P->ctx, P->params, "dtlon").f * DEG_TO_RAD; // rot about +Z = -tlon
        Q->fwd_rot_z[0][0] = cos(theta_z);
        Q->fwd_rot_z[0][1] = -sin(theta_z);
        Q->fwd_rot_z[0][2] = 0;
        Q->fwd_rot_z[1][0] = sin(theta_z);
        Q->fwd_rot_z[1][1] = cos(theta_z);
        Q->fwd_rot_z[1][2] = 0;
        Q->fwd_rot_z[2][0] = 0;
        Q->fwd_rot_z[2][1] = 0;
        Q->fwd_rot_z[2][2] = 1;
        Q->target_pt_rot = true;
    } else {
        for(int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Q->fwd_rot_z[i][j] = i == j ? 1 : 0;
    }
    if (pj_param(P->ctx, P->params, "ttlat").i) {
        double theta_y = M_HALFPI + pj_param(P->ctx, P->params, "dtlat").f * DEG_TO_RAD; // rot about +Z = tlat + pi/2
        Q->fwd_rot_y[0][0] = cos(theta_y);
        Q->fwd_rot_y[0][1] = 0;
        Q->fwd_rot_y[0][2] = sin(theta_y);
        Q->fwd_rot_y[1][0] = 0;
        Q->fwd_rot_y[1][1] = 1;
        Q->fwd_rot_y[1][2] = 0;
        Q->fwd_rot_y[2][0] = -sin(theta_y);
        Q->fwd_rot_y[2][1] = 0;
        Q->fwd_rot_y[2][2] = cos(theta_y);  
        Q->target_pt_rot = true;
    } else {
        for(int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Q->fwd_rot_y[i][j] = i == j ? 1 : 0;
    }

    if (pj_param(P->ctx, P->params, "tstretch_factor").i) {
        Q->stretch_factor = pj_param(P->ctx, P->params, "dstretch_factor").f ;
    } else {
        Q->stretch_factor = 1.0;
    }

    P->inv = gcs_s_inverse;
    P->fwd = gcs_s_forward;
    P->es = 0.;
    P->destructor = gcs_destructor;

    return P;
}
