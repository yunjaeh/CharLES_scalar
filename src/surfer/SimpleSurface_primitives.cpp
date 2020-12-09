#include "SimpleSurface.hpp"
#include "GeomUtils.hpp"
#include "GeodesicSphere.hpp"
#include "SplineStuff.hpp"

int SimpleSurface::addGridForFlaggedTris(const double diam,const double spacing,const double xc[3],const double normal[3],const string& zone_name) {

  // adds a turbulence grid of spheres within (inside and touching) the flagged tris...

  sp_flag.setLength(nsp);

  const double nmag = MAG(normal); assert(nmag > 0.0);
  const double unit_normal[3] = { normal[0]/nmag, normal[1]/nmag, normal[2]/nmag };
  for (int isp = 0; isp < nsp; ++isp) {
    const double dx[3] = DIFF(xsp[isp],xc);
    const double dist = DOT_PRODUCT(dx,unit_normal);
    if (dist > 0.5*diam) {
      sp_flag[isp] = 1;
    }
    else if (dist < -0.5*diam) {
      sp_flag[isp] = -1;
    }
    else {
      sp_flag[isp] = 0;
    }
  }

  int nst_active = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist]) {
      const int sp_sum = sp_flag[spost[ist][0]] + sp_flag[spost[ist][1]] + sp_flag[spost[ist][2]];
      if ((sp_sum != 3)&&(sp_sum != -3)) {
        ++nst_active;
      }
    }
  }

  if (nst_active == 0) {
    cout << "Error: no tris bounding requested grid plane" << endl;
    return -1;
  }

  // use the smallest value in normal to determine a direction not aligned with normal...
  double e1[3];
  if (fabs(normal[0]) < max(fabs(normal[1]),fabs(normal[2]))) {
    // use (1,0,0) cross normal for e1...
    e1[0] = 0.0;
    e1[1] = -normal[2];
    e1[2] = normal[1];
  }
  else if (fabs(normal[1]) < fabs(normal[2])) {
    // use (0,1,0) cross normal for e1...
    e1[0] = normal[2];
    e1[1] = 0.0;
    e1[2] = -normal[0];
  }
  else {
    // use (0,0,1) cross normal for e1...
    e1[0] = -normal[1];
    e1[1] = normal[0];
    e1[2] = 0.0;
  }
  const double e1_mag = MAG(e1);
  FOR_I3 e1[i] /= e1_mag;
  const double e0[3] = CROSS_PRODUCT(unit_normal,e1);

  double bb_min[2] = { 1.0E+20, 1.0E+20 };
  double bb_max[2] = { -1.0E+20, -1.0E+20 };
  int * active_tris = new int[nst_active];
  nst_active = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist]) {
      const int sp_sum = sp_flag[spost[ist][0]] + sp_flag[spost[ist][1]] + sp_flag[spost[ist][2]];
      if ((sp_sum != 3)&&(sp_sum != -3)) {
        active_tris[nst_active++] = ist;
        // also use this active tri to enlarge the bounding box...
        double dist[3];
        FOR_I3 {
          const int isp = spost[ist][i];
          const double dx[3] = DIFF(xsp[isp],xc);
          dist[i] = DOT_PRODUCT(dx,unit_normal);
        }
        FOR_I3 {
          const int isp0 = spost[ist][i];
          const int isp1 = spost[ist][(i+1)%2];
          if (dist[i]*dist[(i+1)%3] < 0.0) {
            double xp[3];
            FOR_J3 xp[j] = (dist[(i+1)%3]*xsp[isp0][j]-dist[i]*xsp[isp1][j])/(dist[(i+1)%3]-dist[i]);
            double x_primed = (xp[0]-xc[0])*e0[0] + (xp[1]-xc[1])*e0[1] + (xp[2]-xc[2])*e0[2];
            double y_primed = (xp[0]-xc[0])*e1[0] + (xp[1]-xc[1])*e1[1] + (xp[2]-xc[2])*e1[2];
            bb_min[0] = min(bb_min[0],x_primed);
            bb_max[0] = max(bb_max[0],x_primed);
            bb_min[1] = min(bb_min[1],y_primed);
            bb_max[1] = max(bb_max[1],y_primed);
          }
          else {
            if (dist[i] == 0.0) {
              double xp[3];
              FOR_J3 xp[j] = xsp[isp0][j];
              double x_primed = (xp[0]-xc[0])*e0[0] + (xp[1]-xc[1])*e0[1] + (xp[2]-xc[2])*e0[2];
              double y_primed = (xp[0]-xc[0])*e1[0] + (xp[1]-xc[1])*e1[1] + (xp[2]-xc[2])*e1[2];
              bb_min[0] = min(bb_min[0],x_primed);
              bb_max[0] = max(bb_max[0],x_primed);
              bb_min[1] = min(bb_min[1],y_primed);
              bb_max[1] = max(bb_max[1],y_primed);
            }
            if (dist[(i+1)%3] == 0.0) {
              double xp[3];
              FOR_J3 xp[j] = xsp[isp1][j];
              double x_primed = (xp[0]-xc[0])*e0[0] + (xp[1]-xc[1])*e0[1] + (xp[2]-xc[2])*e0[2];
              double y_primed = (xp[0]-xc[0])*e1[0] + (xp[1]-xc[1])*e1[1] + (xp[2]-xc[2])*e1[2];
              bb_min[0] = min(bb_min[0],x_primed);
              bb_max[0] = max(bb_max[0],x_primed);
              bb_min[1] = min(bb_min[1],y_primed);
              bb_max[1] = max(bb_max[1],y_primed);
            }
          }
        }
      }
    }
  }

  // here we could introduce an adt to store the tris and find the nearest, but
  // for now, do it brute force...

  // now figure out the grid of points...

  const int i0 = (int)ceil(bb_min[0]/spacing);
  const int i1 = (int)floor(bb_max[0]/spacing);
  const int j0 = (int)ceil(bb_min[1]/spacing);
  const int j1 = (int)floor(bb_max[1]/spacing);

  // check the n*m algorithm...
  if ((i1-i0+1)*(j1-j0+1)*nst_active > 100000) {
    cout << "Warning: this GRID seems very fine: " << (i1-i0+1)*(j1-j0+1) << " spheres, " << nst_active << " active tris (may take a while)." << endl;
  }

  //const double mpoint_factor = 1.0E-3;
  //vector<pair<int,int> > ijVec;
  vector<double> xNoIntersectVec;
  vector<double> xIntersectVec;
  const double diam_tol = 0.05;
  for (int i = i0; i <= i1; ++i) {
    for (int j = j0; j <= j1; ++j) {
      double xp[3];
      FOR_K3 xp[k] = xc[k] + double(i)*spacing*e0[k] + double(j)*spacing*e1[k];
      // for each point, find its closest tri and then determine if the closest
      // tri normal is pointing away from the line of sight. If it is, then
      // we are inside the requested tris...
      int ist_closest = -1;
      double d2_min;
      for (int ist_active = 0; ist_active < nst_active; ++ist_active) {
        const int ist = active_tris[ist_active];
        const double this_d2 = MiscUtils::getPointToTriDist2(xp,xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
        if ((ist_closest == -1)||(this_d2 < d2_min)) {
          ist_closest = ist;
          d2_min = this_d2;
        }
      }
      assert(ist_closest != -1);
      const double tri_normal[3] = TRI_NORMAL_2(xsp[spost[ist_closest][0]],xsp[spost[ist_closest][1]],xsp[spost[ist_closest][2]]);
      const double dx[3] = DIFF(xsp[spost[ist_closest][0]],xp);
      const double dp = DOT_PRODUCT(tri_normal,dx);
      if (dp > 0.0) {
        // dp > 0 means the center of the sphere is inside the geom. It may still intersect the geometry,
        // so treat these cases here...
        if (d2_min > diam*diam/4.0*(1.0+diam_tol)) {
          // this point is completely inside, so some of these points may intersect...
          xNoIntersectVec.push_back(xp[0]);
          xNoIntersectVec.push_back(xp[1]);
          xNoIntersectVec.push_back(xp[2]);
        }
        else if (d2_min > diam*diam/4.0) { // radius stupid!
          // move this point inside so it does not intersect the surface...
          const double delta = sqrt(diam*diam/4.0*(1.0+diam_tol)) - sqrt(d2_min);
          const double mag_n = MAG(tri_normal);
          xNoIntersectVec.push_back(xp[0]-delta*tri_normal[0]/mag_n);
          xNoIntersectVec.push_back(xp[1]-delta*tri_normal[1]/mag_n);
          xNoIntersectVec.push_back(xp[2]-delta*tri_normal[2]/mag_n);
        }
        else {
          // this sphere is touching, but make sure it is touching enough...
          const double delta = max(0.0,sqrt(d2_min) - sqrt(diam*diam/4.0*(1.0-diam_tol)));
          const double mag_n = MAG(tri_normal);
          xIntersectVec.push_back(xp[0]+delta*tri_normal[0]/mag_n);
          xIntersectVec.push_back(xp[1]+delta*tri_normal[1]/mag_n);
          xIntersectVec.push_back(xp[2]+delta*tri_normal[2]/mag_n);
        }
      }
      else if (d2_min < diam*diam/4.0*(1.0-diam_tol)) {
        xIntersectVec.push_back(xp[0]);
        xIntersectVec.push_back(xp[1]);
        xIntersectVec.push_back(xp[2]);
      }
    }
  }

  delete[] active_tris;

  const int n_sphere_edge = 4; // the decimation of the tri...
  int nsp0 = nsp;
  int nst0 = nst;

  const int nsp_inc = GeodesicSphere::getSphereNodeCount(n_sphere_edge);
  const int nst_inc = GeodesicSphere::getSphereTriCount(n_sphere_edge);
  nsp += (xNoIntersectVec.size()+xIntersectVec.size())/3*nsp_inc;
  nst += (xNoIntersectVec.size()+xIntersectVec.size())/3*nst_inc;
  growNspData(nsp,nsp0);
  growNstData(nst,nst0);

  // need a new (or existing) zone for the internal spheres...

  int izone;
  int nzn0 = zoneVec.size();
  for (izone = 0; izone < nzn0; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  // and the bool ones -- put these in the same zone as above...

  /*
  const string zone_name_bool = zone_name+"_bool";
  int izone_bool;
  nzn0 = zoneVec.size();
  for (izone_bool = 0; izone_bool < nzn0; ++izone_bool) {
    if (zoneVec[izone_bool].getName() == zone_name_bool)
      break;
  }
  if (izone_bool == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name_bool));
    ++nsz;
  }
  */
  int izone_bool = izone;

  // NOTE: addSphere assumes memory is available and returns incremented nsp0,nst0...
  for (int ii = 0,limit=xNoIntersectVec.size(); ii < limit; ii += 3) {
    double xp[3];
    xp[0] = xNoIntersectVec[ii];
    xp[1] = xNoIntersectVec[ii+1];
    xp[2] = xNoIntersectVec[ii+2];
    GeodesicSphere::addSphere(xsp+nsp0,spost+nst0,xp,0.5*diam,n_sphere_edge,true);
    for (int ist = nst0; ist < nst0+nst_inc; ++ist) {
      znost[ist] = izone;
      FOR_I3 spost[ist][i] += nsp0;
    }
    nsp0 += nsp_inc;
    nst0 += nst_inc;
  }

  for (int ii = 0,limit=xIntersectVec.size(); ii < limit; ii += 3) {
    double xp[3];
    xp[0] = xIntersectVec[ii];
    xp[1] = xIntersectVec[ii+1];
    xp[2] = xIntersectVec[ii+2];
    GeodesicSphere::addSphere(xsp+nsp0,spost+nst0,xp,0.5*diam,n_sphere_edge,true);
    for (int ist = nst0; ist < nst0+nst_inc; ++ist) {
      znost[ist] = izone_bool;
      FOR_I3 spost[ist][i] += nsp0;
    }
    nsp0 += nsp_inc;
    nst0 += nst_inc;
  }

  assert(nsp0 == nsp);
  assert(nst0 == nst);

  cout << " > added " << xNoIntersectVec.size()+xIntersectVec.size() << " turbulence grid spheres to zone \"" << zone_name << "\"" << endl;

  return 0;

}

int SimpleSurface::addSphere(const double xc0[3],const double r,const int n,const bool flip) {

  // current surface values
  int nsp0 = nsp;
  int nst0 = nst;

  // recall n contains the edge split size: e.g. 1,2,3 or some relatively low int
  nsp += GeodesicSphere::getSphereNodeCount(n);
  nst += GeodesicSphere::getSphereTriCount(n);
  growNspData(nsp,nsp0);
  growNstData(nst,nst0);

  GeodesicSphere::addSphere(xsp+nsp0,spost+nst0,xc0,r,n,flip);

  const string zone_name = "SIMPLE_SPHERE";
  int izone;
  int nzn0 = zoneVec.size();
  for (izone = 0; izone < nzn0; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  for (int ist = nst0; ist < nst; ++ist) {
    znost[ist] = izone;
    FOR_I3 spost[ist][i] += nsp0;
  }

  return 0;

}

int SimpleSurface::addHemisphere(const double xp[3],const double np[3],const double rp,const int ntheta,const bool flip) {

  int nsp0 = nsp;
  int nst0 = nst;

  int nsp_inc,nst_inc;
  GeomUtils::getHemisphereNodeAndTriCount(nsp_inc,nst_inc,ntheta);

  nsp += nsp_inc;
  nst += nst_inc;

  growNspData(nsp,nsp0);
  growNstData(nst,nst0);

  GeomUtils::addHemisphere(xsp+nsp0,spost+nst0,xp,np,rp,ntheta,flip);

  // finally the zone and also offset the spost if we weren't
  // the first tris...
  string zone_name = "HEMISPHERE";
  int izone,nzn;
  for (izone = 0,nzn=zoneVec.size(); izone < nzn; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == int(zoneVec.size())) {
    zoneVec.push_back(SurfaceZone(zone_name));
  }

  for (int ist = nst0; ist < nst; ++ist) {
    znost[ist] = izone;
    FOR_I3 spost[ist][i] += nsp0;
  }

  return 0;

}

int SimpleSurface::addPistonFF(const double dr0,const double dz0,const double dr1,const double dz1) {

  assert(0); // march 2019
  // probably use lifted surface and cylinder trimming...

  ensureTeost();

  // put the edges in a vec...

  int * sp_flag = new int[nsp];
  FOR_ISP sp_flag[isp] = -1;

  vector<pair<int,int> > edgeVec;
  double xcc[3] = { 0.0, 0.0, 0.0 };
  double normal[3] = { 0.0, 0.0, 0.0 };
  double d2min[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
  double d2max[3] = { 0.0, 0.0, 0.0 };
  double wgt_sum = 0.0;
  FOR_IST {
    FOR_I3 {
      if (isEdgeOpen(ist,i)) {
        edgeVec.push_back(pair<int,int>(ist,i));
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        sp_flag[isp0] = 0;
        sp_flag[isp1] = 0;
        const double wgt = DIST(xsp[isp0],xsp[isp1]);
        FOR_I3 xcc[i] += wgt*(xsp[isp0][i]+xsp[isp1][i]);
        wgt_sum += wgt;
        const double cp[3] = CROSS_PRODUCT(xsp[isp0],xsp[isp1]);
        FOR_I3 normal[i] += cp[i];
        // also, assuming the piston is centered in either x,y,z, compute the min
        // and max radius^2...
        double xmid[3] = { 0.0, 0.0, 0.0 };
        FOR_I3 {
          xmid[i] += 0.5*(xsp[isp0][i]+xsp[isp1][i]);
          d2max[i] = max(d2max[i],xsp[isp0][(i+1)%3]*xsp[isp0][(i+1)%3]+xsp[isp0][(i+2)%3]*xsp[isp0][(i+2)%3]);
          d2max[i] = max(d2max[i],xsp[isp1][(i+1)%3]*xsp[isp1][(i+1)%3]+xsp[isp1][(i+2)%3]*xsp[isp1][(i+2)%3]);
        }
        FOR_I3 {
          d2min[i] = min(d2min[i],xmid[(i+1)%3]*xmid[(i+1)%3]+xmid[(i+2)%3]*xmid[(i+2)%3]);
        }
      }
    }
  }
  FOR_I3 xcc[i] /= 2.0*wgt_sum;
  double mag = MAG(normal);
  FOR_I3 normal[i] /= mag;

  cout << " > xcc, normal: " << COUT_VEC(xcc) << " " << COUT_VEC(normal) << endl;

  int nsp_edge = 0;
  FOR_ISP {
    if (sp_flag[isp] == 0)
      sp_flag[isp] = nsp_edge++;
    else
      assert(sp_flag[isp] == -1);
  }

  // if this is a loop, the nsp_edge size should be the same as the edgeVec.size()...
  assert(nsp_edge == int(edgeVec.size()));

  int id = -1;
  if (fabs(normal[0]) > 10.0*max(fabs(normal[1]),fabs(normal[2]))) {
    // this is an x-aligned piston...
    assert(0);
  }
  else if (fabs(normal[1]) >10.0*max(fabs(normal[0]),fabs(normal[2]))) {
    // this is an y-aligned piston...
    assert(0);
  }
  else if (fabs(normal[2]) >10.0*max(fabs(normal[0]),fabs(normal[1]))) {
    // this is an z-aligned piston...
    cout << " > assuming z-alignment" << endl;
    id = 2;
  }
  else {
    // no obvious alignment direction...
    assert(0);
  }

  cout << " > min/max radius of open edges: " << sqrt(d2min[id]) << " " << sqrt(d2max[id]) << endl;

  // the number of new tris will be edgeVec.size*5, and new nodes is 2*nsp_edge+1 (for the center)...
  const int nst0 = nst;
  nst += edgeVec.size()*5;
  growNstData(nst,nst0);

  const int nsp0 = nsp;
  nsp += nsp_edge*2+1;
  growNspData(nsp,nsp0);

  const string zone_name = "PISTON_FF";
  int izone;
  int nzn0 = zoneVec.size();
  for (izone = 0; izone < nzn0; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  // z-only for now...
  assert(id == 2);

  // nodes first...
  for (int isp = 0; isp < nsp0; ++isp) {
    if (sp_flag[isp] >= 0) {
      // figure out our "r"...
      const double rmag = sqrt(xsp[isp][0]*xsp[isp][0] + xsp[isp][1]*xsp[isp][1]);
      // first ring of points...
      xsp[nsp0+sp_flag[isp]][0] = xsp[isp][0] + dr0*xsp[isp][0]/rmag;
      xsp[nsp0+sp_flag[isp]][1] = xsp[isp][1] + dr0*xsp[isp][1]/rmag;
      xsp[nsp0+sp_flag[isp]][2] = xsp[isp][2] + dz0;
      // second ring of points...
      xsp[nsp0+nsp_edge+sp_flag[isp]][0] = xsp[isp][0] + dr1*xsp[isp][0]/rmag;
      xsp[nsp0+nsp_edge+sp_flag[isp]][1] = xsp[isp][1] + dr1*xsp[isp][1]/rmag;
      xsp[nsp0+nsp_edge+sp_flag[isp]][2] = xsp[isp][2] + dz1;
    }
  }
  assert(nsp0+nsp_edge*2 == nsp-1);
  xsp[nsp0+nsp_edge*2][0] = xcc[0];
  xsp[nsp0+nsp_edge*2][1] = xcc[1];
  xsp[nsp0+nsp_edge*2][2] = xcc[2]+dz1;

  // then new tris...
  for (int ied = 0,ned=edgeVec.size(); ied < ned; ++ied) {
    const int ist = edgeVec[ied].first;
    const int i = edgeVec[ied].second;
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    assert(sp_flag[isp0] >= 0);
    assert(sp_flag[isp1] >= 0);
    // first pair of tris...
    spost[nst0+ied*5  ][0] = isp1;
    spost[nst0+ied*5  ][1] = isp0;
    spost[nst0+ied*5  ][2] = nsp0+sp_flag[isp1];
    spost[nst0+ied*5+1][0] = isp0;
    spost[nst0+ied*5+1][1] = nsp0+sp_flag[isp0];
    spost[nst0+ied*5+1][2] = nsp0+sp_flag[isp1];
    // second pair of tris...
    spost[nst0+ied*5+2][0] = nsp0+sp_flag[isp1];
    spost[nst0+ied*5+2][1] = nsp0+sp_flag[isp0];
    spost[nst0+ied*5+2][2] = nsp0+nsp_edge+sp_flag[isp1];
    spost[nst0+ied*5+3][0] = nsp0+sp_flag[isp0];
    spost[nst0+ied*5+3][1] = nsp0+nsp_edge+sp_flag[isp0];
    spost[nst0+ied*5+3][2] = nsp0+nsp_edge+sp_flag[isp1];
    // third tri...
    spost[nst0+ied*5+4][0] = nsp0+nsp_edge+sp_flag[isp1];
    spost[nst0+ied*5+4][1] = nsp0+nsp_edge+sp_flag[isp0];
    spost[nst0+ied*5+4][2] = nsp0+nsp_edge*2;
    // znost...
    znost[nst0+ied*5  ] = izone;
    znost[nst0+ied*5+1] = izone;
    znost[nst0+ied*5+2] = izone;
    znost[nst0+ied*5+3] = izone;
    znost[nst0+ied*5+4] = izone;
  }

  delete[] sp_flag;

  return 0;

}

int SimpleSurface::addOffset(const double dn) {

  ensureTeost();

  // put the edges in a vec...

  vector<pair<int,int> > edgeVec;
  double xcc[3] = { 0.0, 0.0, 0.0 };
  double normal[3] = { 0.0, 0.0, 0.0 };
  double d2min[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
  double d2max[3] = { 0.0, 0.0, 0.0 };
  double wgt_sum = 0.0;
  FOR_IST {
    FOR_I3 {
      if (isEdgeOpen(ist,i)) {
        edgeVec.push_back(pair<int,int>(ist,i));
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        const double wgt = DIST(xsp[isp0],xsp[isp1]);
        FOR_I3 xcc[i] += wgt*(xsp[isp0][i]+xsp[isp1][i]);
        wgt_sum += wgt;
        const double cp[3] = CROSS_PRODUCT(xsp[isp0],xsp[isp1]);
        FOR_I3 normal[i] += cp[i];
        // also, assuming the piston is centered in either x,y,z, compute the min
        // and max radius^2...
        double xmid[3] = { 0.0, 0.0, 0.0 };
        FOR_I3 {
          xmid[i] += 0.5*(xsp[isp0][i]+xsp[isp1][i]);
          d2max[i] = max(d2max[i],xsp[isp0][(i+1)%3]*xsp[isp0][(i+1)%3]+xsp[isp0][(i+2)%3]*xsp[isp0][(i+2)%3]);
          d2max[i] = max(d2max[i],xsp[isp1][(i+1)%3]*xsp[isp1][(i+1)%3]+xsp[isp1][(i+2)%3]*xsp[isp1][(i+2)%3]);
        }
        FOR_I3 {
          d2min[i] = min(d2min[i],xmid[(i+1)%3]*xmid[(i+1)%3]+xmid[(i+2)%3]*xmid[(i+2)%3]);
        }
      }
    }
  }
  FOR_I3 xcc[i] /= 2.0*wgt_sum;
  double mag = MAG(normal);
  FOR_I3 normal[i] /= mag;

  cout << " > xcc, normal: " << COUT_VEC(xcc) << " " << COUT_VEC(normal) << endl;

  // if this is a loop, the nsp_edge size should be the same as the edgeVec.size()...
  //assert(nsp_edge == edgeVec.size());

  // the number of new tris will be edgeVec.size*5, and new nodes is 2*nsp_edge+1 (for the center)...
  const int nst0 = nst;
  nst += nst0 + edgeVec.size()*2;
  growNstData(nst,nst0);

  const int nsp0 = nsp;
  nsp += nsp0;
  growNspData(nsp,nsp0);

  const string zone_name = "OFFSET";
  int izone;
  int nzn0 = zoneVec.size();
  for (izone = 0; izone < nzn0; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  // build the normals at the old nodes...
  double (*n_sp)[3] = new double[nsp0][3];
  for (int isp = 0; isp < nsp0; ++isp)
    FOR_I3 n_sp[isp][i] = 0.0;

  for (int ist = 0; ist < nst0; ++ist) {
    const double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
    FOR_I3 {
      const int isp = spost[ist][i];
      FOR_J3 n_sp[isp][j] += n[j];
    }
  }

  // normalize...
  for (int isp = 0; isp < nsp0; ++isp) {
    const double mag = MAG(n_sp[isp]);
    assert(mag > 0.0);
    FOR_I3 n_sp[isp][i] /= mag;
  }

  // place new nodes at old + dn*n_sp...
  for (int isp = 0; isp < nsp0; ++isp) {
    const int isp_new = isp+nsp0;
    FOR_I3 xsp[isp_new][i] = xsp[isp][i] + n_sp[isp][i]*dn;
  }

  // then new tris: tri copies first: note orientation must be flipped...
  for (int ist = 0; ist < nst0; ++ist) {
    const int ist_new = ist+nst0;
    spost[ist_new][0] = spost[ist][0]+nsp0;
    spost[ist_new][1] = spost[ist][2]+nsp0;
    spost[ist_new][2] = spost[ist][1]+nsp0;
    znost[ist_new] = izone;
  }

  // then any tris from open edges...
  for (int ied = 0,ned=edgeVec.size(); ied < ned; ++ied) {
    const int ist = edgeVec[ied].first;
    const int i = edgeVec[ied].second;
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    // pair of tris...
    // 1...
    spost[nst0*2+ied*2][0] = isp0;
    spost[nst0*2+ied*2][1] = isp0+nsp0;
    spost[nst0*2+ied*2][2] = isp1;
    znost[nst0*2+ied*2]    = izone;
    // 2...
    spost[nst0*2+ied*2+1][0] = isp1;
    spost[nst0*2+ied*2+1][1] = isp0+nsp0;
    spost[nst0*2+ied*2+1][2] = isp1+nsp0;
    znost[nst0*2+ied*2+1]    = izone;
  }

  delete[] n_sp;
  return 0;

}

int SimpleSurface::constructPlane(int& iarg,Param * param) {
  int ierr = -1;
  const double pt[3] = {param->getDouble(iarg),param->getDouble(iarg+1),param->getDouble(iarg+2)};
  iarg += 3;
  const double norm[3] = {param->getDouble(iarg),param->getDouble(iarg+1),param->getDouble(iarg+2)};
  iarg += 3;

  int nx = 1;
  int ny = 1;
  double width,height = -1.0;

  while (iarg < param->size()) {
    const string tok = param->getString(iarg++);
    if (tok == "WIDTH") {
      width = param->getDouble(iarg++);
    }
    else if (tok == "HEIGHT") {
      height = param->getDouble(iarg++);
    }
    else if (tok == "NX") {
      nx = param->getInt(iarg++);
    }
    else if (tok == "NY") {
      ny = param->getInt(iarg++);
    }
    else {
      CWARN("unrecognized parameter \"" << tok << "\"; skipping");
    }
  }

  if (width<=0.0 || height <= 0.0) {
    CWARN("a positive WIDTH and HEIGHT must be specified");
    ierr = -1;
  }
  else ierr = addPlane(pt,norm,width,height,nx,ny);

  return ierr;
}

// initialize a plane surface by specifying the center point, normal axis, width, and height
int SimpleSurface::addPlane(const double xp[3],const double np[3],const double width,const double height,const int _nx, const int _ny) {

  // make sure enough points to discretize edges
  const int nx = max(_nx,1);
  const int ny = max(_ny,1);

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += (nx+1)*(ny+1);  // number of nodes is +1 number of facets in each dimension
  growNspData(nsp,ss_nsp0);

  double up[3] = {0.0,0.0,0.0};
  if (fabs(np[0]) < min(fabs(np[1]),fabs(np[2]))  || (np[0] == 0.0 && np[1] == 0.0)) up[0] = 1.0;
  else if (fabs(np[1]) < min(fabs(np[0]),fabs(np[2]))) up[1] = 1.0;
  else up[2] = 1.0;

  double e0[3] = CROSS_PRODUCT(up,np);
  NORMALIZE(e0);
  double e1[3] = CROSS_PRODUCT(np,e0);
  NORMALIZE(e1);

  const double delta_w = fabs(width)/nx;
  const double delta_h = fabs(height)/ny;

  double top_left[3];
  FOR_I3 top_left[i] = xp[i] - 0.5*fabs(height)*e0[i] - 0.5*fabs(width)*e1[i];

  int count = 0;
  for (int row=0, nrow=ny+1; row<nrow; ++row) {
    for (int col=0,ncol=nx+1; col<ncol; ++col) {
      ++count;
      FOR_I3 xsp[ss_nsp0 + (row*ncol) + col][i] = top_left[i] + row*delta_h*e0[i] + col*delta_w*e1[i];
    }
  }
  assert(count == (nsp-ss_nsp0));

  const string zone_name = "SIMPLE_PLANE";
  int new_zone;
  int nzn0 = zoneVec.size();
  for (new_zone = 0; new_zone < nzn0; ++new_zone) {
    if (zoneVec[new_zone].getName() == zone_name)
      break;
  }
  if (new_zone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  nst += 2*nx*ny;  // number of new tris
  growNstData(nst,ss_nst0);

  count = 0;
  for (int row=0; row<ny; ++row) {
    for (int col=0; col<nx; ++col) {
      // below indices are relative to new nodes; need to index appropriately
      const int ul = row*(nx+1) + col;
      const int ur = ul+1;
      const int ll = ul+nx+1;
      const int lr = ll+1;

      // add 2 tris for this square
      spost[ss_nst0 + 2*(row*nx + col)][0] = ss_nsp0 + ul;
      spost[ss_nst0 + 2*(row*nx + col)][1] = ss_nsp0 + ll;
      spost[ss_nst0 + 2*(row*nx + col)][2] = ss_nsp0 + lr;

      spost[ss_nst0 + 2*(row*nx + col)+1][0] = ss_nsp0 + ul;
      spost[ss_nst0 + 2*(row*nx + col)+1][1] = ss_nsp0 + lr;
      spost[ss_nst0 + 2*(row*nx + col)+1][2] = ss_nsp0 + ur;

      count+=2;
    }
  }
  assert(count == nst-ss_nst0);

  for (int ist=ss_nst0; ist<nst; ++ist) {
    znost[ist] = new_zone;
  }

  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;
}

int SimpleSurface::addCircle(const double xp[3],const double np[3],const double r,const int n) {
  assert(n>=3);
  assert(MAG(np)>0.0);

  double e1[3];
  double e2[3];
  MiscUtils::getBestE1E2FromE0(e1,e2,np);

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += n+1;  // number of nodes is n+1
  growNspData(nsp,ss_nsp0);

  FOR_I3 xsp[ss_nsp0][i] = xp[i];
  MiscUtils::createCirclePts(xsp,ss_nsp0+1,xp,np,r,n);

  // populate tri information
  const int new_zone = zoneVec.size();
  zoneVec.push_back(SurfaceZone("circle"));

  nst += n;
  growNstData(nst,ss_nst0);

  facetCircleToPoint(spost,znost,ss_nsp0+0,ss_nsp0+1,ss_nst0+0,new_zone,n,false);

  return 0;
}

int SimpleSurface::addAnnulus(const double xp[3],const double np[3],const double r0,const double r1,const int n) {
  assert(n>=3);
  assert(MAG(np)>0.0);

  double e1[3];
  double e2[3];
  MiscUtils::getBestE1E2FromE0(e1,e2,np);

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += 2*n;
  growNspData(nsp,ss_nsp0);

  MiscUtils::createCirclePts(xsp,ss_nsp0,xp,np,r0,n);
  MiscUtils::createCirclePts(xsp,ss_nsp0+n,xp,np,r1,n);

  // populate tri information
  const int new_zone = zoneVec.size();
  zoneVec.push_back(SurfaceZone("annulus"));

  nst += 2*n;
  growNstData(nst,ss_nst0);

  facetCircleToCircle(spost,znost,ss_nsp0+0,ss_nsp0+n,ss_nst0+0,new_zone,n,true,true);

  return 0;
}

// initialize a box surface from 2 corner points...
int SimpleSurface::addBox(const double x0[3],const double x1[3],const bool flip) {

  // make a box consisting of 6 sides, and 12 tris...
  // hex: standard node numbering is...
  //
  //        7------6
  //       /      /|
  //      4------5 |
  //      | 3    | 2
  //      |      |/
  //      0------1
  //

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += 8;
  growNspData(nsp,ss_nsp0);

  xsp[ss_nsp0+0][0] = x0[0]; xsp[ss_nsp0+0][1] = x0[1]; xsp[ss_nsp0+0][2] = x0[2];
  xsp[ss_nsp0+1][0] = x1[0]; xsp[ss_nsp0+1][1] = x0[1]; xsp[ss_nsp0+1][2] = x0[2];
  xsp[ss_nsp0+2][0] = x1[0]; xsp[ss_nsp0+2][1] = x1[1]; xsp[ss_nsp0+2][2] = x0[2];
  xsp[ss_nsp0+3][0] = x0[0]; xsp[ss_nsp0+3][1] = x1[1]; xsp[ss_nsp0+3][2] = x0[2];
  xsp[ss_nsp0+4][0] = x0[0]; xsp[ss_nsp0+4][1] = x0[1]; xsp[ss_nsp0+4][2] = x1[2];
  xsp[ss_nsp0+5][0] = x1[0]; xsp[ss_nsp0+5][1] = x0[1]; xsp[ss_nsp0+5][2] = x1[2];
  xsp[ss_nsp0+6][0] = x1[0]; xsp[ss_nsp0+6][1] = x1[1]; xsp[ss_nsp0+6][2] = x1[2];
  xsp[ss_nsp0+7][0] = x0[0]; xsp[ss_nsp0+7][1] = x1[1]; xsp[ss_nsp0+7][2] = x1[2];

  // do this a little more robustly in the future in case there are matching zone
  // names already there...
  const int zone_x0 = zoneVec.size();
  const int zone_x1 = zoneVec.size()+1;
  const int zone_y0 = zoneVec.size()+2;
  const int zone_y1 = zoneVec.size()+3;
  const int zone_z0 = zoneVec.size()+4;
  const int zone_z1 = zoneVec.size()+5;
  zoneVec.push_back(SurfaceZone("x0"));
  zoneVec.push_back(SurfaceZone("x1"));
  zoneVec.push_back(SurfaceZone("y0"));
  zoneVec.push_back(SurfaceZone("y1"));
  zoneVec.push_back(SurfaceZone("z0"));
  zoneVec.push_back(SurfaceZone("z1"));
  nsz += 6;

  nst += 12;
  growNstData(nst,ss_nst0);

  // x0...
  spost[ss_nst0+0][0] = 0; spost[ss_nst0+0][1] = 4; spost[ss_nst0+0][2] = 7; znost[ss_nst0+0] = zone_x0;
  spost[ss_nst0+1][0] = 0; spost[ss_nst0+1][1] = 7; spost[ss_nst0+1][2] = 3; znost[ss_nst0+1] = zone_x0;

  // x1...
  spost[ss_nst0+2][0] = 1; spost[ss_nst0+2][1] = 6; spost[ss_nst0+2][2] = 5; znost[ss_nst0+2] = zone_x1;
  spost[ss_nst0+3][0] = 1; spost[ss_nst0+3][1] = 2; spost[ss_nst0+3][2] = 6; znost[ss_nst0+3] = zone_x1;

  // y0...
  spost[ss_nst0+4][0] = 0; spost[ss_nst0+4][1] = 1; spost[ss_nst0+4][2] = 5; znost[ss_nst0+4] = zone_y0;
  spost[ss_nst0+5][0] = 0; spost[ss_nst0+5][1] = 5; spost[ss_nst0+5][2] = 4; znost[ss_nst0+5] = zone_y0;

  // y1...
  spost[ss_nst0+6][0] = 2; spost[ss_nst0+6][1] = 3; spost[ss_nst0+6][2] = 6; znost[ss_nst0+6] = zone_y1;
  spost[ss_nst0+7][0] = 3; spost[ss_nst0+7][1] = 7; spost[ss_nst0+7][2] = 6; znost[ss_nst0+7] = zone_y1;

  // z0...
  spost[ss_nst0+8][0] = 0; spost[ss_nst0+8][1] = 3; spost[ss_nst0+8][2] = 2; znost[ss_nst0+8] = zone_z0;
  spost[ss_nst0+9][0] = 0; spost[ss_nst0+9][1] = 2; spost[ss_nst0+9][2] = 1; znost[ss_nst0+9] = zone_z0;

  // z1...
  spost[ss_nst0+10][0] = 4; spost[ss_nst0+10][1] = 5; spost[ss_nst0+10][2] = 6; znost[ss_nst0+10] = zone_z1;
  spost[ss_nst0+11][0] = 4; spost[ss_nst0+11][1] = 6; spost[ss_nst0+11][2] = 7; znost[ss_nst0+11] = zone_z1;

  if (ss_nsp0 > 0) {
    // if adding surface, properly offset node indices, zones
    for (int ist = ss_nst0; ist < nst; ++ist) {
      FOR_I3 spost[ist][i] += ss_nsp0;
    }
  }

  // if flip requested...
  if (flip) {
    for (int ist = ss_nst0; ist < nst; ++ist) {
      const int tmp = spost[ist][1];
      spost[ist][1] = spost[ist][2];
      spost[ist][2] = tmp;
    }
  }

  return 0;

}

int SimpleSurface::addDisk(const double x0[3],const double x1[3],const double r,const int n,const bool b_flip) {
  const int ss_nst0 = nst;
  const int ierr = addTcone(x0,x1,r,r,n);
  if ((ierr == 0) && b_flip) {
    for (int ist = ss_nst0; ist < nst; ++ist) {
      const int tmp = spost[ist][1];
      spost[ist][1] = spost[ist][2];
      spost[ist][2] = tmp;
    }
  }
  return ierr;
}

int SimpleSurface::addTcone(const double x0[3],const double x1[3],const double rad0,const double rad1,const int n) {

  // n points around base...
  const double dx[3] = DIFF(x1,x0);
  const double dx_mag = MAG(dx);
  if (dx_mag <= 0.0) {
    CERR("cannot create TCONE with overlapping endpoints");
  }

  // absolute value ensures no singularity in the middle of the tcone as well as positive volume
  const double r0 = fabs(rad0);
  const double r1 = fabs(rad1);

  const double SMALL = 1E-12;
  if (r0 < SMALL || r1 < SMALL) {
    CERR("cannot create a TCONE with either radius = 0");
  }

  double axis[3];
  FOR_I3 axis[i] = x1[i] - x0[i];
  COUT2(" TCONE axis = " << COUT_VEC( axis ) );

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += 2 + 2*n;
  growNspData(nsp,ss_nsp0);

  // first and last node are cap centers
  FOR_I3 {
    xsp[ss_nsp0+0][i] = x0[i];
    xsp[nsp - 1][i] = x1[i];
  }

  // create nodes around each cap, cap0 first [1:n],
  // then cap1 [n+1:2n]
  // TODO: modify to GeomUtils...
  MiscUtils::createCirclePts(xsp,ss_nsp0+1,x0,axis,r0,n);
  MiscUtils::createCirclePts(xsp,ss_nsp0+n+1,x1,axis,r1,n);

  const string zone_name = "tcone";
  /*
  int new_zone;
  int nzn0 = zoneVec.size();
  for (new_zone = 0; new_zone < nzn0; ++new_zone) {
    if (zoneVec[new_zone].getName() == zone_name)
      break;
  }
  if (new_zone == nzn0) {
  */
  const int new_zone = zoneVec.size();
  zoneVec.push_back(SurfaceZone(zone_name));
  ++nsz;
  //}

  nst += n*4;
  growNstData(nst,ss_nst0);

  // cap0, cap1
  // TODO: and here...
  facetCircleToPoint(spost,znost,ss_nsp0+0,ss_nsp0+1,ss_nst0+0,new_zone,n,false);
  facetCircleToPoint(spost,znost,nsp-1,ss_nsp0+n+1,ss_nst0+3*n,new_zone,n,true);

  // wall
  facetCircleToCircle(spost,znost,ss_nsp0+1,ss_nsp0+n+1,ss_nst0+n,new_zone,n,true,true);

  return 0;
}

int SimpleSurface::addAnnularTcone(const double x0[3],const double x1[3],const double rad00,const double rad01,const double rad10,const double rad11,const int n) {

  // n points around base...
  const double dx[3] = DIFF(x1,x0);
  const double dx_mag = MAG(dx);
  if (dx_mag <= 0.0) {
    CERR("cannot create ANNULAR_TCONE with overlapping endpoints");
  }

  // absolute value ensures no singularity in the middle of the tcone as well as positive volume
  const double r00 = fabs(rad00);
  const double r01 = fabs(rad01);
  const double r10 = fabs(rad10);
  const double r11 = fabs(rad11);

  const double SMALL = 1E-12;
  if (r00 < SMALL || r01 < SMALL || r10 < SMALL || r11 < SMALL) {
    CERR("cannot create an ANNUULAR_TCONE with any radius = 0");
  }

  double axis[3];
  FOR_I3 axis[i] = x1[i] - x0[i];
  COUT2(" ANNULAR_TCONE axis = " << COUT_VEC( axis ) );

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  nsp += 4*n;
  growNspData(nsp,ss_nsp0);


  // create nodes around each cap
  MiscUtils::createCirclePts(xsp,ss_nsp0+0*n,x0,axis,r00,n);
  MiscUtils::createCirclePts(xsp,ss_nsp0+1*n,x0,axis,r01,n);

  MiscUtils::createCirclePts(xsp,ss_nsp0+2*n,x1,axis,r10,n);
  MiscUtils::createCirclePts(xsp,ss_nsp0+3*n,x1,axis,r11,n);

  const string zone_name = "annular_tcone";
  int new_zone;
  int nzn0 = zoneVec.size();
  for (new_zone = 0; new_zone < nzn0; ++new_zone) {
    if (zoneVec[new_zone].getName() == zone_name)
      break;
  }
  if (new_zone == nzn0) {
    zoneVec.push_back(SurfaceZone(zone_name));
    ++nsz;
  }

  nst += n*8;
  growNstData(nst,ss_nst0);

  // cap0
  facetCircleToCircle(spost,znost,ss_nsp0+0*n,ss_nsp0+1*n,ss_nst0+0*n,new_zone,n,true,true);
  // r1
  facetCircleToCircle(spost,znost,ss_nsp0+1*n,ss_nsp0+3*n,ss_nst0+2*n,new_zone,n,true,true);
  // r0
  facetCircleToCircle(spost,znost,ss_nsp0+2*n,ss_nsp0+0*n,ss_nst0+4*n,new_zone,n,true,true);
  // cap1
  facetCircleToCircle(spost,znost,ss_nsp0+3*n,ss_nsp0+2*n,ss_nst0+6*n,new_zone,n,true,true);

  return 0;
}

int SimpleSurface::addRevolve(const double rAxis[3],const double x0[3],double (*ir_vals)[2],const int n_vals,const int _n_theta, bool capSurf, const double maxAR) {
  COUT1("SimpleSurface::initRevolve()");

  if (n_vals < 2) {
    CWARN("not enough points were specified in the profile; skipping");
    return -1;
  }
  if (_n_theta < 3) {
    CWARN("N_THETA of \"" << _n_theta << "\" is too small; defaulting to 3");
  }
  const int n_theta = max(3,_n_theta);

  double axis[3];
  FOR_I3 {
    axis[i] = rAxis[i];
  }
  NORMALIZE(axis); // make a unit vector
  COUT2(" > revolve axis = " << COUT_VEC( axis ) );

  vector< pair<double,double> > irVec;
  for (int v=0; v < n_vals; ++v) {
    irVec.push_back(pair<double,double> (ir_vals[v][0],ir_vals[v][1]));
  }

  for (int ixy = 1, limit=irVec.size()-1; ixy < limit; ++ixy) {
    if (irVec[ixy].second <= 0.0) {
      CWARN("cannot support profile points on the axis itself unless the starting or ending point; skipping");
      return -1;
    }
  }

  bool specifiedCap0 = false;
  bool specifiedCap1 = false;
  if (irVec[0].second == 0.0) {
    specifiedCap0 = true;
  } else {
    assert(irVec[0].second > 0.0);
  }

  if (irVec[irVec.size()-1].second == 0.0) {
    specifiedCap1 = true;
  } else {
    assert(irVec[irVec.size()-1].second > 0.0);
  }

  // if specified, add points within irVec to maintain certain aspect ratio
  // of facets. Since we use a fixed "n", we only look at adjustments introduced
  // by adding points in (I,R) space.
  if (maxAR != -1.0) {
    assert(maxAR >= 1.0);
    COUT2(" > aspect ratio limit: " << maxAR);

    if (capSurf) {
      // add point(s) to close the surface
      if (!specifiedCap0) {
        irVec.insert( irVec.begin(), pair<double,double>(irVec.front().first,0.0) );
        COUT2(" > adding starting point to surface at: (" << irVec.front().first << "," << irVec.front().second << ")");
      }

      if (!specifiedCap1) {
        irVec.push_back( pair<double,double>(irVec.back().first,0.0) );
        COUT2(" > adding closing point to surface at: (" << irVec.back().first << "," << irVec.back().second << ")");
      }
    }

    vector< pair<double,double> > irVec_new;
    irVec_new.push_back(irVec[0]);  //start from first point in IR profile

    int pos = 0;
    for (int i=1,limit=irVec.size(); i < limit; i++) {
      // COUT2("segment[" << i << "]");
      const double dI = irVec[i].first - irVec_new[pos].first;
      const double dR = irVec[i].second - irVec_new[pos].second;
      const double L0 = sqrt(dI*dI + dR*dR);
      double L1;
      if (i != limit-1) {
        L1 = 2*irVec[i].second*sin(M_PI/double(n_theta));
      }
      else {
        L1 = 2*irVec_new[pos].second*sin(M_PI/double(n_theta));
      }
      const double L2 = sqrt(L0*L0 + L1*L1);

      // compute aspect ratio based on ratio of circumradius/(2*inradius)
      // if ratio below criteria, add point(s) to profile
      const double s = 0.5*(L0+L1+L2);
      const double AR = L0*L1*L2 / (8.0*(s-L0)*(s-L1)*(s-L2));

      // COUT2(" > lengths L0,L1,L2,s: " << L0 << " " << L1 << " " << L2 << " " << s);
      // COUT2(" > computed AR: " << AR);

      if ((AR > maxAR) && (max(L0,L2) > L1)) {
        const int newPoints = (int)floor(AR/maxAR);
        // COUT2("    > inserting " << newPoints << " points");
        const double deltaI = dI / double(newPoints+1);
        const double deltaR = dR / double(newPoints+1);
        for (int p=0; p < newPoints; p++) {
          const double newI = irVec_new[pos].first + (p+1)*deltaI;
          const double newR = irVec_new[pos].second + (p+1)*deltaR;
          irVec_new.push_back(pair<double,double>(newI,newR));
        }
        pos += newPoints;
      }
      // add specified point
      irVec_new.push_back(irVec[i]);
      ++pos;
    }
    irVec.swap(irVec_new);

    if (capSurf) {
      // remove endpoints from irVec - used only for AR limiting
      if (!specifiedCap0) irVec.erase(irVec.begin());
      if (!specifiedCap1) irVec.pop_back();
    }
  }

  double x1[3];
  if (capSurf) {
    FOR_I3 x1[i] = x0[i] + irVec.back().first*axis[i];
  }

  if (specifiedCap0) {
    // remove endpoints from irVecs
    irVec.erase(irVec.begin());
  }

  if (specifiedCap1) {
    FOR_I3 x1[i] = x0[i] + irVec.back().first*axis[i];
    irVec.pop_back();
  }
  COUT2("x0:" << COUT_VEC(x0));
  COUT2("x1:" << COUT_VEC(x1));

  for (int i=0,limit=irVec.size(); i < limit; i++) {
    COUT2("[" << i << "]: " << irVec[i].first << " " << irVec[i].second);
  }
  COUT2(" > number of IR points prescribed: " << irVec.size());


  // Add new points
  const int ss_nsp0 = nsp;

  nsp += irVec.size()*n_theta;
  if (capSurf || specifiedCap0) nsp += 1;  // account for starting point
  if (capSurf || specifiedCap1) nsp += 1;  // account for ending point

  growNspData(nsp,ss_nsp0);

  int isp = ss_nsp0;
  if (capSurf || specifiedCap0) {
    FOR_I3 xsp[isp][i] = x0[i];
    isp++;  // cap0 endpoint
  }

  // place all irVec points
  for (int ixy=0, limit=irVec.size(); ixy < limit; ++ixy) {
    double center[3];
    FOR_I3 {
      center[i] = x0[i] + irVec[ixy].first*axis[i];
    }
    MiscUtils::createCirclePts(xsp,isp,center,axis,irVec[ixy].second,n_theta,false); // do NOT stagger the points
    isp += n_theta;
  }

  if (capSurf || specifiedCap1) {
    FOR_I3 xsp[isp][i] = x1[i];
    isp++;
  }

  assert(isp == nsp);

  // Add new tris
  const int ss_nst0 = nst;

  nst += 2*n_theta*(irVec.size()-1);
  if (capSurf || specifiedCap0) nst += n_theta;  // accounts for starting cap facets
  if (capSurf || specifiedCap1) nst += n_theta;  // accounts for ending cap facets

  growNstData(nst,ss_nst0);
  const int ss_nzn0 = zoneVec.size();

  COUT1(" > zones:");
  zoneVec.push_back(SurfaceZone("wall")); COUT1("    > \"wall\"");
  int cap0_idx=0;  // default write to wall zone
  int cap1_idx=0;  // default write to wall zone
  if (capSurf || specifiedCap0) {
    zoneVec.push_back(SurfaceZone("cap0")); COUT1("    > \"cap0\"");
    cap0_idx++;
    cap1_idx++;
  }
  if (capSurf || specifiedCap1) {
    zoneVec.push_back(SurfaceZone("cap1")); COUT1("    > \"cap1\"");
    cap1_idx++;
  }

  // place all surface tris
  int ist = ss_nst0;

  if (capSurf || specifiedCap0) {
    // cap0
    facetCircleToPoint(spost,znost,ss_nsp0,ss_nsp0+1,ist,cap0_idx,n_theta,false);
    ist += n_theta;
  }

  for (int ixy=0,limit=irVec.size()-1; ixy < limit; ++ixy) {
    // revolve
    int c0  = ss_nsp0 + ixy*n_theta;
    int c1  = ss_nsp0 + ixy*n_theta + n_theta;
    if (capSurf || specifiedCap0) {
      c0 += 1;
      c1 += 1;
    }
    facetCircleToCircle(spost,znost,c0,c1,ist,0,n_theta,true,true);
    ist += 2*n_theta;
  }

  if (capSurf || specifiedCap1) {
    // cap1
    facetCircleToPoint(spost,znost,nsp-1,nsp-n_theta-1,ist,cap1_idx,n_theta,true);
    ist += n_theta;
  }

  assert(ist == nst);

  if (ss_nst0 > 0) {
    // if adding surface, properly offset zone indices
    for (int ist = ss_nst0; ist < nst; ++ist) {
      // FOR_I3 spost[ist][i] += ss_nsp0;
      znost[ist] += ss_nzn0;
    }
  }

  COUT1(" > addRevolve done");
  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;  // success
}

int SimpleSurface::addExtrudedProfile(const double x0[3], const double x1[3], double (*xyz_vals)[3], const int n_vals, const int _n_span, const bool closed) {

  COUT1("SimpleSurface::addExtrudedProfile()");
  COUT1(" > assumes profile is closed loop; if not this will produce errors currently");
  if (n_vals < 2) {
    CWARN("not enough points were specified in the profile; skipping");
    return -1;
  }

  int n_span = (_n_span+1);  // number of profile stations total, including start and end
  if (_n_span < 1) {
    CWARN("invalid specification of N_SPAN; defaulting to N_SPAN = 1");
    n_span = 2;
  }

  const double axis[3] = DIFF(x1,x0);
  if (MAG(axis) <= 0.0) {
    CWARN("extrude direction is not properly defined; skipping");
    return -1;
  }
  COUT2(" > extrusion direction: " << COUT_VEC(axis));

  // Add new points
  const int ss_nsp0 = nsp;
  nsp += n_vals*n_span;

  growNspData(nsp,ss_nsp0);

  // fill in points
  const double denom = 1.0 / (n_span-1);
  for (int ispan=0; ispan < n_span; ++ispan) {
    const double frac = ispan*denom;
    // double xc[3];
    // FOR_I3 xc[i] = x0[i] + frac*axis[i];
    for (int ii=0; ii<n_vals; ++ii) {
      FOR_I3 xsp[ss_nsp0+(ispan*n_vals)+ii][i] = x0[i] + xyz_vals[ii][i] + frac*axis[i];
    }
    // MiscUtils::create3DFrom2DPtAndNormal(xsp,ispan*n_vals,xc,axis,xy_vals,n_vals);
  }

  // create all tris
  const int ss_nst0 = nst;

  int _n_vals = n_vals;
  if (!closed) --_n_vals;  // last panel (between end and start node) doesn't get generated
  nst += 2*(n_span-1)*_n_vals;

  growNstData(nst,ss_nst0);
  const int ss_nzn0 = zoneVec.size();
  COUT1(" > zones:");
  zoneVec.push_back(SurfaceZone("extrusion")); COUT1("    > \"extrusion\"");

  int ist = ss_nst0;
  for (int ispan=0, limit = n_span-1; ispan < limit; ++ispan) {
    // start indices for circles-of-points nodes
    const int c0 = ss_nsp0 + ispan*n_vals;
    const int c1 = ss_nsp0 + (ispan+1)*n_vals;
    facetCircleToCircle(spost,znost,c0,c1,ist,ss_nzn0,n_vals,closed,true);
    ist += 2*_n_vals;
  }

  assert(ist == nst);

  COUT1(" > addExtrude done");
  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;  // success

}

// create the facets between a circle of points and a single point

void SimpleSurface::facetCircleToPoint(int (* spost)[3], int * znost, const int indexPt, const int indexCircle, const int st0, const int zoneId, const int n, const bool flip) {
  // flip dictates node ordering (normal direction convention)
  int pt1_inc = 1;
  int pt2_inc = 0;
  if (flip) {
    pt1_inc = 0;
    pt2_inc = 1;
  }

  for (int i=0, limit=n; i < limit; ++i) {
    spost[st0+i][0] = indexPt;
    spost[st0+i][1] = indexCircle + (i+pt1_inc)%n;
    spost[st0+i][2] = indexCircle + (i+pt2_inc)%n;
    znost[st0+i] = zoneId;
  }
}


// create the facets between two circles of points, assuming they have the same n

void SimpleSurface::facetCircleToCircle(int (* spost)[3], int * znost, const int indexC0, const int indexC1, const int st0, const int zoneId, const int n, const bool closed, const bool flip) {
  // flip dictates node ordering (normal direction convention)
  int pt1_inc = 1;
  int pt2_inc = 0;
  if (flip) {
    pt1_inc = 0;
    pt2_inc = 1;
  }

  int limit = n;
  if (!closed) --limit;
  for (int i=0; i < limit; ++i) {
    spost[st0+(2*i)][0] = indexC0 + (i+pt1_inc)%n;
    spost[st0+(2*i)][1] = indexC0 + (i+pt2_inc)%n;
    spost[st0+(2*i)][2] = indexC1 + (i+pt1_inc)%n;
    znost[st0+(2*i)] = zoneId;

    spost[st0+(2*i + 1)][0] = indexC0 + (i+pt2_inc)%n;;
    spost[st0+(2*i + 1)][1] = indexC1 + (i+pt2_inc)%n;
    spost[st0+(2*i + 1)][2] = indexC1 + (i+pt1_inc)%n;
    znost[st0+(2*i + 1)] = zoneId;
  }
}

int SimpleSurface::addLofted(const vector<string>& profileVec,const int nr) {

  using namespace SplineStuff;

  const int npr = profileVec.size();
  assert(npr >= 2);
  vector<double> * dVec = new vector<double>[npr];
  int np = -1; // number of points in all profiles (should be the same)...
  for (int ipr = 0; ipr < npr; ++ipr) {
    //cout << " > reading BLADEGEN profile: " << profileVec[ipr] << endl;
    GeomUtils::readXYZ(dVec[ipr],profileVec[ipr]);
    assert(dVec[ipr].size()%3 == 0);
    const int this_np = dVec[ipr].size()/3;
    //cout << " > profile has np: " << this_np << endl;
    if (np == -1) {
      np = this_np;
    }
    else if (np != this_np) {
      cout << "Error: profiles with different point counts not supported" << endl;
      return -1;
    }
  }

  CubicSpline * cspline = new CubicSpline[np];
  double (*xpr)[3] = new double[npr][3];
  for (int ip = 0; ip < np; ++ip) {
    for (int ipr = 0; ipr < npr; ++ipr) {
      FOR_I3 xpr[ipr][i] = dVec[ipr][ip*3+i];
    }
    cspline[ip].init(xpr,npr);
  }
  delete[] xpr;

  // the number of tris is np*nr*2, and the number of nodes is np*(nr+1)...

  const int nsp = np*(nr+1);
  double (*xsp)[3] = new double[nsp][3];
  int isp = 0;
  for (int ip = 0; ip < np; ++ip) {
    for (int ir = 0; ir <= nr; ++ir) {
      cspline[ip].getX(xsp[isp],double(ir)/double(nr));
      ++isp;
    }
  }
  assert(isp == nsp);

  const int nst = np*nr*2;
  int (*spost)[3] = new int[nst][3];
  int ist = 0;
  for (int ip0 = 0; ip0 < np; ++ip0) {
    const int ip1 = (ip0+1)%np;
    for (int ir = 0; ir < nr; ++ir) {
      spost[ist][0] = ip0*(nr+1)+ir;
      spost[ist][1] = ip1*(nr+1)+ir;
      spost[ist][2] = ip0*(nr+1)+ir+1;
      spost[ist+1][0] = ip1*(nr+1)+ir;
      spost[ist+1][1] = ip1*(nr+1)+ir+1;
      spost[ist+1][2] = ip0*(nr+1)+ir+1;
      ist += 2;
    }
  }
  assert(ist == nst);

  // for now, just output an sbin...

  GeomUtils::writeTecplot("nasa37.dat",spost,nst,xsp);
  GeomUtils::writeSbin("nasa37.sbin",spost,nst,xsp);

  cout << " > looks good: about to throw" << endl;
  getchar();
  throw(0);

}

int SimpleSurface::addExtrudeSpline(vector<double>& dVec,const double x0[3],const double x1[3],const int nspan) {

  using namespace SplineStuff;

  // the dVec contains x0 and x1 locations...

  assert(dVec.size()%2 == 0);
  const int np = dVec.size()/2;

  // splines are 3D, so use z = 0...

  double (*xp)[3] = new double[np][3];
  for (int ip = 0; ip < np; ++ip) {
    xp[ip][0] = dVec[ip*2  ];
    xp[ip][1] = dVec[ip*2+1];
    xp[ip][2] = 0.0; // 2d spline would be nice
  }
  CubicSpline cspline;
  cspline.init(xp,np,true);
  delete[] xp;
  const double length = cspline.getLength();

  cout << " > length: " << length << endl;

  // take a look...
  {
    FILE * fp = fopen("spline.dat","w");
    const int n = 2000;
    for (int i = 0; i <= n; ++i) {
      double xr[3];
      cspline.getX(xr,double(i)/double(n));
      fprintf(fp,"%18.15le %18.15le\n",xr[0],xr[1]);
    }
    fclose(fp);
  }

  // to determine nx, we use aspect ratio considerations...

  const double dx[3] = DIFF(x1,x0);
  const double dx_mag = MAG(dx);
  assert(dx_mag > 0.0);
  double e1[3],e2[3];
  MiscUtils::getBestE1E2FromE0(e1,e2,dx);

  const int nx = max(1,int(length/dx_mag*double(nspan)));
  cout << " > nspan: " << nspan << ", nx (computed from isotropic considerations): " << nx << endl;

  int nsp = nx*(nspan+1);
  int nst = 2*nspan*nx;

  cout << " > nsp: " << nsp << ", nst: " << nst << endl;

  double (*xsp)[3] = new double[nsp][3];
  int isp = 0;
  for (int i = 0; i < nx; ++i) {
    double x[3];
    cspline.getX(x,double(i)/double(nx));
    for (int j = 0; j <= nspan; ++j) {
      FOR_K3 xsp[isp][k] = x0[k] + double(j)/double(nspan)*dx[k] + x[0]*e1[k] + x[1]*e2[k];
      ++isp;
    }
  }
  assert(isp == nsp);

  int ist = 0;
  int (*spost)[3] = new int[nst][3];
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nspan; ++j) {
      // this patch has 4 points:
      const int isp0 = i*(nspan+1)+j;
      const int isp1 = ((i+1)%nx)*(nspan+1)+j;
      const int isp2 = ((i+1)%nx)*(nspan+1)+(j+1);
      const int isp3 = i*(nspan+1)+(j+1);
      spost[ist][0] = isp0;
      spost[ist][1] = isp1;
      spost[ist][2] = isp2;
      spost[ist+1][0] = isp0;
      spost[ist+1][1] = isp2;
      spost[ist+1][2] = isp3;
      ist += 2;
    }
  }
  assert(ist == nst);

  //GeomUtils::writePtsTecplot("pts.dat",xsp,nsp);
  GeomUtils::writeSbin("extrude_spline.sbin",spost,nst,xsp);
  GeomUtils::writeTecplot("extrude_spline.dat",spost,nst,xsp);

  delete[] spost;
  delete[] xsp;

  cout << "EXTRUDE_SPLINE just wrote extrude_spline.sbin -- finish this off some day" << endl;

  return -1;

}

int SimpleSurface::addRevolveSpline(vector<double>& dVec,const double axis[3],const double origin[3],const int ntheta) {

  using namespace SplineStuff;

  // the dVec contains axial distance, radial distance...

  assert(dVec.size()%2 == 0);
  const int np = dVec.size()/2;

  // splines are 3D, so use z = 0...

  double r_avg = 0.0;
  double (*xp)[3] = new double[np][3];
  for (int ip = 0; ip < np; ++ip) {
    xp[ip][0] = dVec[ip*2  ];
    xp[ip][1] = dVec[ip*2+1];
    r_avg += xp[ip][1];
    xp[ip][2] = 0.0; // 2d spline would be nice
  }
  CubicSpline cspline;
  cspline.init(xp,np);
  delete[] xp;
  r_avg /= double(np);
  const double length = cspline.getLength();

  cout << " > r_avg: " << r_avg << " length: " << length << endl;

  // take a look...
  {
    FILE * fp = fopen("spline.dat","w");
    const int n = 2000;
    for (int i = 0; i <= n; ++i) {
      double xr[3];
      cspline.getX(xr,double(i)/double(n));
      fprintf(fp,"%18.15le %18.15le\n",xr[0],xr[1]);
    }
    fclose(fp);
  }

  const double axis_mag = MAG(axis);
  assert(axis_mag > 0.0);
  double e1[3],e2[3];
  MiscUtils::getBestE1E2FromE0(e1,e2,axis);

  // to determine nx, we use aspect ratio considerations:

  const int nx = max(1,int(length/(2.0*M_PI*r_avg/double(ntheta))));
  cout << " > ntheta: " << ntheta << ", nx (computed from isotropic considerations): " << nx << endl;

  int nsp = (nx+1)*ntheta;
  int nst = 2*ntheta*nx;

  cout << " > nsp: " << nsp << ", nst: " << nst << endl;

  double (*xsp)[3] = new double[nsp][3];
  int isp = 0;
  for (int i = 0; i <= nx; ++i) {
    double xr[3];
    cspline.getX(xr,double(i)/double(nx));
    for (int j = 0; j < ntheta; ++j) {
      const double theta = M_PI*2.0*double(j)/double(ntheta);
      FOR_K3 xsp[isp][k] = origin[k] + xr[0]*axis[k] + xr[1]*(cos(theta)*e1[k] + sin(theta)*e2[k]);
      ++isp;
    }
  }
  assert(isp == nsp);

  int ist = 0;
  int (*spost)[3] = new int[nst][3];
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ntheta; ++j) {
      // this patch has 4 points:
      const int isp0 = i*ntheta+j;
      const int isp1 = (i+1)*ntheta+j;
      const int isp2 = (i+1)*ntheta+((j+1)%ntheta);
      const int isp3 = i*ntheta+((j+1)%ntheta);
      spost[ist][0] = isp0;
      spost[ist][1] = isp1;
      spost[ist][2] = isp2;
      spost[ist+1][0] = isp0;
      spost[ist+1][1] = isp2;
      spost[ist+1][2] = isp3;
      ist += 2;
    }
  }
  assert(ist == nst);

  //GeomUtils::writePtsTecplot("pts.dat",xsp,nsp);
  GeomUtils::writeSbin("revolve_spline.sbin",spost,nst,xsp);
  GeomUtils::writeTecplot("revolve_spline.dat",spost,nst,xsp);

  delete[] spost;
  delete[] xsp;

  cout << "REVOLVE_SPLINE just wrote revolve_spline.sbin -- finish this off some day" << endl;
  return -1;

}
