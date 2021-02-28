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

static PJ_XY gcs_s_forward (PJ_LP lp, PJ *P) {           /* Spheroidal, forward */
    PJ_XY xy = {0.0,0.0};
    struct pj_opaque *Q = static_cast<struct pj_opaque*>(P->opaque);
    double temp;
    enum Face i;
    PJ_LP lp_i;

    int s_end = Q->selected_face ? 1 : 6;  // if a face is selected it is first in MRU list, t.f., stop search for s>=1

    for(int s = 0; s < s_end; ++s) {
        i = static_cast<enum Face>(Q->mru_face[s]);
        lp_i = {lp.lam - Q->face_false_origin[i].lam, lp.phi};
        xy = Q->face_proj[i]->fwd(lp_i, Q->face_proj[i]);
        // correct face if proj errno is 0, and -1 <= x <= 1, and -1 <= y <= 1
        if (proj_errno_reset(Q->face_proj[i]) == 0 && xy.x >= -1 && xy.x <= 1 && xy.y >= -1 && xy.y <= 1) {
            // Update MRU list
            for(; s >= 1; --s) {
                Q->mru_face[s] = Q->mru_face[s-1];
            }
            Q->mru_face[0] = i;
            break;
        }
    }

    if (xy.x > 1 || xy.x < -1 || xy.y > 1 || xy.y < -1) {
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

    // Initialize MRU list
    for (int i = X_P; i <= Z_N; ++i) {
        Q->mru_face[i] = Q->specific_face ? Q->selected_face : static_cast<enum Face>(i);
    }

    P->inv = gcs_s_inverse;
    P->fwd = gcs_s_forward;
    P->es = 0.;
    P->destructor = gcs_destructor;

    return P;
}
