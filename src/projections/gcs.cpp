#define PJ_LIB__

#include <errno.h>
#include <math.h>

#include "proj.h"
#include "proj_internal.h"
#include <math.h>

PROJ_HEAD(gcs, "GnominicCubedSphere") "\n\tAzi, Sph";

#define EPS10  1.e-10

PJ_XY decode_tile_xy(char c, const char** error_msg);
void get_tile_rotations(PJ_XY* tile_pos, int* ccw_rotations, const char** error_msg);

namespace { // anonymous namespace
constexpr double bn_trans_mat_xyz_to_uv[6][2][3] = { // xyz->uv transformation matrix for the base net (bn); a base net is one that requires no rotations
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
} // anonymous namespace

namespace { // anonymous namespace
struct pj_opaque {
    int tile_rotation_mat[6][2][2];
    int selected_face;
    bool explicit_face;
    int face_search_list[6];
    PJ_XY tile_false_origin[6];
    PJ_LP face_false_origin[6];
    PJ* face_proj[6];
    char face_proj_def[6][256];
};
} // anonymous namespace

PJ_XY decode_tile_xy(char c, const char** error_msg) {
    PJ_XY xy = {0.0,0.0};
    int v = toupper(c);
    *error_msg = nullptr;

    // decode character 
    if (v >= '0' && v <= '9') {
        v -= '0';
    } else if (v >= 'A' && v <= 'H') {
        v -= 'A';
        v += 10;
    } else if (v >= 'J' && v <= 'K') {
        v -= 'J';
        v += 18;
    } else {
        *error_msg = "Invalid character in +cube_net (valid characters are 0-9,A-H,J-K case-insensitive)";
    }

    // XY offset
    if (v >= 0 && v <= 15) {
        xy.x = static_cast<double>(v % 4);
        xy.y = -static_cast<double>(v / 4);
    } else if (v == 16) {
        xy.x = 0.;
        xy.y = -4.;
    } else if (v == 17) {
        xy.x = 1.;
        xy.y = -4.;
    } else if (v == 18) {
        xy.x = 4.;
        xy.y = 0.;
    } else if (v == 19) {
        xy.x = 4.;
        xy.y = -1.;
    } else {
        *error_msg = "Internal error decoding tile offset";
    }

    // multiply by 2
    xy.x *= 2.;
    xy.y *= 2.;
    return xy;
}

void get_tile_rotations(PJ_XY* tile_pos, int* ccw_rotations, const char** error_msg) {
    error_msg = nullptr;
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
    int i, s;
    PJ_LP lp_i;
    for(s = 0; s < 6; ++s) {
        i = Q->face_search_list[s];
        lp_i = {lp.lam - Q->face_false_origin[i].lam, lp.phi};
        xy = Q->face_proj[i]->fwd(lp_i, Q->face_proj[i]);
        if (proj_errno_reset(Q->face_proj[i]) == 0 && xy.x >= -1 && xy.x <= 1 && xy.y >= -1 && xy.y <= 1) {
            break;
        }
    }

    while(s > 0) {
        Q->face_search_list[s] = Q->face_search_list[s-1];
        s--;
    }
    Q->face_search_list[0] = i;


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

    for(i = 0; i < 6; ++i) {
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

// need destructor

PJ *PROJECTION(gcs) {
    struct pj_opaque *Q = static_cast<struct pj_opaque*>(calloc (1, sizeof (struct pj_opaque)));
    PJ_XY px_tile_pos; // position of +x tile; other tile position are relative to +x
    int cube_index[6] = {0, 3, 1, 4, 2, 5}; // {+x, -x, +y, -y, +z, -z}, default value is GEOS-5 indexing
    const char* err_msg = nullptr;
    if (nullptr==Q)
        return pj_default_destructor (P, PROJ_ERR_OTHER /*ENOMEM*/);
    P->opaque = Q;

    // get index of faces
    if (pj_param(P->ctx, P->params, "tcube_index").i) {
        const char* arg_cube_index = pj_param(P->ctx, P->params, "scube_index").s;

        if (strlen(arg_cube_index) != 6) {
            proj_log_error(P, _("+cube_index must be exactly 6 characters"));
            return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }

        for(int i = 0; i < 6; ++i) {
            cube_index[i] = static_cast<int>(arg_cube_index[i] - '0');
            if (cube_index[i] < 0 || cube_index[i] > 6) {
                proj_log_error(P, _("Invalid index in +cube_index. Indices must be 0-6."));
                return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
            }
        }
    }

    // get explicit face if specified
    Q->explicit_face = pj_param(P->ctx, P->params, "tcube_face").i;
    if (Q->explicit_face) {
        int cube_face = pj_param(P->ctx, P->params, "icube_face").i;
        Q->selected_face = -1;
        for(int i = 0; i < 6; ++i) {
            if (cube_face == cube_index[i]) {
                Q->selected_face = i;
            }
        }
        if (Q->selected_face == -1) {
            proj_log_error(P, _("Invalid +cube_face (value doesn't occur in +cube_index)."));
            return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }
    }

    // get cube net
    if (pj_param(P->ctx, P->params, "tcube_net").i) {
        const char* cube_net = pj_param(P->ctx, P->params, "scube_net").s;

        if (strlen(cube_net) != 6) {
            proj_log_error(P, _("+cube_net must be exactly 6 characters"));
            return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }

        px_tile_pos = decode_tile_xy(cube_net[0], &err_msg);
        if (err_msg != nullptr) {
            proj_log_error(P, _(err_msg));
            return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }
        Q->tile_false_origin[0] = {0.0, 0.0};

        for (int i = 1; i < 6; ++i) {
            Q->tile_false_origin[i] = decode_tile_xy(cube_net[i], &err_msg);
            Q->tile_false_origin[i].x -= px_tile_pos.x;
            Q->tile_false_origin[i].y -= px_tile_pos.y;
            if (err_msg != nullptr) {
                proj_log_error(P, _(err_msg));
                return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
            }
        }
    } else {
        /* Default cube index */
        Q->tile_false_origin[0] = {0.0, 0.0};
        px_tile_pos = decode_tile_xy('8', &err_msg); // +x
        Q->tile_false_origin[1] = decode_tile_xy('6', &err_msg); // -x
        Q->tile_false_origin[2] = decode_tile_xy('9', &err_msg); // +y
        Q->tile_false_origin[3] = decode_tile_xy('2', &err_msg); // -y
        Q->tile_false_origin[4] = decode_tile_xy('5', &err_msg); // +z
        Q->tile_false_origin[5] = decode_tile_xy('3', &err_msg); // +z

        for(int i = 1; i < 6; ++i) {
            Q->tile_false_origin[i].x -= px_tile_pos.x;
            Q->tile_false_origin[i].y -= px_tile_pos.y;
        }
    }

    // get tile rotations
    int ccw_rotations[6];
    get_tile_rotations(Q->tile_false_origin, ccw_rotations, &err_msg);
    if (err_msg != nullptr) {
        proj_log_error(P, _(err_msg));
        return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
    }

    // tile rotation matrices
    for (int i = 0; i < 6; ++i) {
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

    // projections for each face
    Q->face_false_origin[0] = {0., 0.};
    Q->face_false_origin[1] = {M_PI, 0.};
    Q->face_false_origin[2] = {M_HALFPI, 0.};
    Q->face_false_origin[3] = {-M_HALFPI, 0.};
    Q->face_false_origin[4] = {0., M_HALFPI};
    Q->face_false_origin[5] = {0., -M_HALFPI};
    for (int i = 0; i < 6; ++i) {
        sprintf(Q->face_proj_def[i], "+proj=gnom +lat_0=%f +lon_0=%f", Q->face_false_origin[i].phi*180./M_PI, Q->face_false_origin[i].lam*180./M_PI);
        Q->face_proj[i] = proj_create(P->ctx, Q->face_proj_def[i]);
        if (Q->face_proj[i] == nullptr) {
            return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
        }
    }

    for (int i = 0; i < 6; ++i) {
        Q->face_search_list[i] = i;
    }

    P->inv = gcs_s_inverse;
    P->fwd = gcs_s_forward;
    P->es = 0.;

    // printf("gcs projection\n");
    // printf("\t-- false origins x: %f, %f, %f, %f, %f, %f\n", 
    //     Q->tile_false_origin[0].x,
    //     Q->tile_false_origin[1].x,
    //     Q->tile_false_origin[2].x,
    //     Q->tile_false_origin[3].x,
    //     Q->tile_false_origin[4].x,
    //     Q->tile_false_origin[5].x
    // );
    // printf("\t-- false origins y: %f, %f, %f, %f, %f, %f\n", 
    //     Q->tile_false_origin[0].y,
    //     Q->tile_false_origin[1].y,
    //     Q->tile_false_origin[2].y,
    //     Q->tile_false_origin[3].y,
    //     Q->tile_false_origin[4].y,
    //     Q->tile_false_origin[5].y
    // );
    // printf("\t-- rotations: %i, %i, %i, %i, %i, %i\n", 
    //     ccw_rotations[0],
    //     ccw_rotations[1],
    //     ccw_rotations[2],
    //     ccw_rotations[3],
    //     ccw_rotations[4],
    //     ccw_rotations[5]
    // );

    return P;
}
