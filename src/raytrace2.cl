__kernel void checkintersection(__global float* nodes,__global float* rays,
__global float* ts, int elem_size)
{
    int _i = blockIdx.x;
    // int _j = get_global_id(1);
    double ax,ay,az,bx,by,bz,e01x,e01y,e01z,e02x,e02y,e02z,nl,px,py,pz,d,d_1,
    tx,ty,tz,u,v,w,qx,qy,qz,modr,nnx,nny,nnz,rdnx,rdny,rdnz, D, Dy, Dz, t;

    double min_dis = 1e20;
    
    rox = rays[_i*6];
    roy = rays[_i*6+1];
    roz = rays[_i*6+2];
        
    rdx = rays[_i*6+3];
    rdy = rays[_i*6+4];
    rdz = rays[_i*6+5];
    
    for(unsigned int _j = 0; _j < (*elem_size)*19; _j+=19)
    {

        ax = nodes[_j];
        ay = nodes[_j+1];
        az = nodes[_j+2];

        bx = nodes[_j+3];
        by = nodes[_j+4];
        bz = nodes[_j+5];

        cx = nodes[_j+6];
        cy = nodes[_j+7];
        cz = nodes[_j+8];

        nx = nodes[_j + 9];
        ny = nodes[_j + 10];
        nz = nodes[_j + 11];

        e01x = nodes[_j + 12];
        e01y = nodes[_j + 13];
        e01z = nodes[_j + 14];
                
        e02x = nodes[_j + 15];
        e02y = nodes[_j + 16];
        e02z = nodes[_j + 17];

        nl = nodes[_j + 18];

        D = - nl*(nx * rdx + ny*rdy + nz*rdz);  // |-d e1 e2| = -n.d

        if(fabs(D) < 1e-6)
        {
            continue;
        }

        tx = rox - ax;
        ty = roy - ay;
        tz = roz - az;
        // double px, py, pz;
        // cross(rdx, rdy, rdz, e02x, e02y, e02z, px, py, pz);

        Dy = rdx*(tz*e02y - ty*e02z) + tx*(rdy*e02z - e02y*rdz) + e02x*(ty*rdz - rdy*tz);
        u = Dy / D;

        if(u < 0 || u > 1)
        {
            continue;
        }

        Dz = rdx*(e01z*ty - e01y*tz) + e01x*(tz*rdy - ty*rdz) + tx*(e01y*rdz - rdy*e01z);

        v = Dz / D;

        if(v < 0 || v + u > 1)
        {
            continue;
        }

        t = (tx*nx + ty*ny + tz*nz)*nl/D;

        if(t < 0)
        {
            continue;
        }
        if(t < min_dis)
        {
            min_dis = t;
        }
    }
    if(min_dis < 1e19)
    {
        ts[_i] = min_dis;
    }
}