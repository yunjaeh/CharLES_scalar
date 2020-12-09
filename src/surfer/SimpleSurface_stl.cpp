#include "SimpleSurface.hpp"
#include "DoubleVertex.hpp"

#define FAST_STL_BUFFER_SIZE 65535
#define FAST_STL_TOKEN_SIZE 65535

char * getNextStlToken(FILE * fp,int& pos,int& max_pos,char * buf) {

  char * token = NULL;
  while (1) {

    if (pos >= max_pos) {

      assert(pos == max_pos);
      if (token) {
        // if the token is active, but not completed, then we need to shift the token part of
        // the current buf (i.e. the end) to the start of the buf, and read the next part in...
        //cout << "pos >= max_pos: " << pos << " " << max_pos << " token[0]: " << token[0] << " max_pos-token+buf: " << max_pos-int(token-buf) << endl;
        pos = max_pos - int(token-buf);
        if (token != buf) {
          memmove(buf,token,pos);
          token = buf; // reset to start...
        }
      }
      else {
        pos = 0;
      }

      max_pos = pos + fread(buf+pos, sizeof(char), FAST_STL_BUFFER_SIZE, fp);
      if (max_pos == pos) {
        buf[pos] = '\0';
        return token;
      }

    }

    //cout << "getNextStlToken pos: " << pos << " c: \"" << buf[pos] << "\"" << endl;
    const char c = buf[pos++];

    if ( (c == ' ') || (c == '\t') || (c == '\n') || (c == 13) ) {
      // c is a whitespace character - this represents either the termination
      // of a token, or some space at the start of the next token...
      if (token) break;
    }
    else if (!token) {
      // any other character is the start of the token...
      token = buf+pos-1;
    }

  }

  // terminate string and return...
  buf[pos-1] = '\0';
  return token;

}

#define STL_HEADER_SIZE 84
#define STL_SIZEOF_FACET 50
#define STL_MIN_FILE_SIZE 284
#define STL_LABEL_SIZE 80

int SimpleSurface::addStl(const string& filename) {

  FILE * fp = NULL;
  const int file_err = MiscUtils::openFile(&fp,filename,"rb");
  if (file_err != 0) return file_err;

  // determine if the file is binary or ascii...
  // skip over binary header and test
  unsigned char chtest[2*STL_HEADER_SIZE];
  fseek(fp, STL_HEADER_SIZE, SEEK_SET);
  int pos = fread(chtest,sizeof(chtest),1,fp);
  bool is_binary = false;
  for(int i = 0, end=sizeof(chtest); i < end; ++i) {
    if(chtest[i] > 127) {
      is_binary = true;
      break;
    }
  }

  // above check can fail for short files, so also check ascii first token
  {
    fseek(fp,0,SEEK_SET);
    int _pos = 0;
    int _max_pos = 0;
    char * buf_tmp = new char[FAST_STL_BUFFER_SIZE+FAST_STL_TOKEN_SIZE];
    char * token = getNextStlToken(fp,_pos,_max_pos,buf_tmp);
    if (strcmp(token,"solid") == 0) is_binary = false;
    DELETE(buf_tmp);
  }

  vector<DoubleVertex> vertexVec;
  vector<int> edgeVec;
  //vector<pair<string,int> > newZoneAndLastIvVec;
  vector<pair<int,int> > newIzoneAndLastIvVec;
  map<const string,int> zoneNameMap;


  if (is_binary) {
    COUT1(" > filetype: binary");

    // for a binary file, we get the number of facets from the file size...

    fseek(fp, 0, SEEK_END);
    const int8 file_size = ftell(fp);
    //if ( ((file_size - STL_HEADER_SIZE)%STL_SIZEOF_FACET != 0) || (file_size < STL_MIN_FILE_SIZE) )
    //We should allow binary stls with almost no facets...
    if ((file_size - STL_HEADER_SIZE) % STL_SIZEOF_FACET != 0) {
      // wrong file size...
      cerr << "file_size = " << file_size << endl;
      cerr << "STL_HEADER_SIZE = " << STL_HEADER_SIZE << endl;
      cerr << "STL_SIZEOF_FACET = " << STL_SIZEOF_FACET << endl;
      cerr << "STL_MIN_FILE_SIZE = " << STL_MIN_FILE_SIZE << endl;
      cerr << "The file has the wrong size " << endl;
      throw(-1);
    }

    const int nfa = (file_size - STL_HEADER_SIZE)/STL_SIZEOF_FACET;
    COUT1(" > filetype: binary");
    COUT1(" > number of facets: " << nfa);

    // go to the end of the header...

    fseek(fp,STL_HEADER_SIZE,SEEK_SET);

    // and start reading the chunks...

    map<const pair<DoubleVertex,DoubleVertex>,int> edgeMap;
    vertexVec.resize(nfa*3);
    edgeVec.resize(nfa*3);

    int nv = 0;
    for (int ifa = 0; ifa < nfa; ++ifa) {

      fread(chtest,sizeof(unsigned char),STL_SIZEOF_FACET,fp);

      // the first 3 floats are the normal, which we can discard...

      const int iv0 = nv++;
      const int iv1 = nv++;
      const int iv2 = nv++;

      //Can no longer copy directly from chtest (float) to
      //vertexVec (double) without converting values.
      //Maybe there is a more efficient way...?
      float readBuf[3];
      memcpy(readBuf,chtest+sizeof(float)*3,sizeof(float)*3);
      vertexVec[iv0].x = readBuf[0]; //float to double
      vertexVec[iv0].y = readBuf[1]; //float to double
      vertexVec[iv0].z = readBuf[2]; //float to double
      memcpy(readBuf,chtest+sizeof(float)*6,sizeof(float)*3);
      vertexVec[iv1].x = readBuf[0]; //float to double
      vertexVec[iv1].y = readBuf[1]; //float to double
      vertexVec[iv1].z = readBuf[2]; //float to double
      memcpy(readBuf,chtest+sizeof(float)*9,sizeof(float)*3);
      vertexVec[iv2].x = readBuf[0]; //float to double
      vertexVec[iv2].y = readBuf[1]; //float to double
      vertexVec[iv2].z = readBuf[2]; //float to double


      // use the flag?...
      //uint2 flag;
      //memcpy(&flag,chtest+sizeof(float)*12,sizeof(uint2));
      //cout << "flag: " << flag << endl;

      // skip collapsed tris
      if ((vertexVec[iv0] == vertexVec[iv0+1])||(vertexVec[iv0] == vertexVec[iv0+2])||(vertexVec[iv0+1] == vertexVec[iv0+2])) {
        // rewind nv...
        nv -= 3;
      }
      else {
        // look for edge pairs...
        for (int i = 0; i < 3; ++i) {
          map<const pair<DoubleVertex,DoubleVertex>,int>::iterator iter = edgeMap.find(pair<DoubleVertex,DoubleVertex>(vertexVec[iv0+i],vertexVec[iv0+(i+1)%3]));

          if (iter != edgeMap.end()) {
            assert(edgeVec[iter->second] == -1);
            edgeVec[iter->second] = iv0+i;
            edgeVec[iv0+i] = iter->second;
            // found a match...
            edgeMap.erase(iter);
          }
          else {
            // not found, so add the reverse edge to edgeMap...
            edgeMap[pair<DoubleVertex,DoubleVertex>(vertexVec[iv0+(i+1)%3],vertexVec[iv0+i])] = iv0+i;
            edgeVec[iv0+i] = -1;
          }
        }
      }

    }

    fclose(fp);

    // we may have skipped some colapsed tris...

    if (nv < nfa*3) {
      assert(nv%3 == 0);
      vertexVec.resize(nv);
      edgeVec.resize(nv);
    }

    // for the newIzoneAndLastIvVec, just put everything in a single zone...
    size_t filename_pos = filename.find_last_of("/\\");
    if (filename_pos!=string::npos) zoneNameMap[filename.substr(filename_pos+1)] = 0;
    else zoneNameMap[filename] = 0;

    newIzoneAndLastIvVec.push_back(pair<int,int>(0,nv));

  }
  else {

    // ------------------------
    // ascii version...
    // ------------------------
    COUT1(" > filetype: ascii");
    // rewind to the beginning of the file...
    fseek(fp,0,SEEK_SET);

    pos = 0;
    int max_pos = 0;
    char * buf = new char[FAST_STL_BUFFER_SIZE+FAST_STL_TOKEN_SIZE];
    //bool got_solid = false;
    //string solid_name;
    map<const pair<DoubleVertex,DoubleVertex>,int> edgeMap;

    int * lastIvPtr = NULL;

    while (char * token = getNextStlToken(fp,pos,max_pos,buf)) {
      //cout << "just got token: \"" << token << "\"." << endl;
      if (strcmp(token,"solid") == 0) {
        const string solid_name = getNextStlToken(fp,pos,max_pos,buf);
        map<const string,int>::iterator iter = zoneNameMap.find(solid_name);
        int izone;
        if (iter == zoneNameMap.end()) {
          cout << " > solid \"" << solid_name << "\"" << endl;
          izone = zoneNameMap.size();
          zoneNameMap[solid_name] = izone;
        }
        else {
          izone = iter->second;
        }
        newIzoneAndLastIvVec.push_back(pair<int,int>(izone,vertexVec.size()));
        lastIvPtr = &(newIzoneAndLastIvVec.back().second);
      }
      else if (strcmp(token,"endsolid") == 0) {
        assert(vertexVec.size()%3 == 0);
        //cout << vertexVec.size()/3 << " tris OK" << endl;
      }
      else if (strcmp(token,"facet") == 0) {
        assert(vertexVec.size()%3 == 0);
      }
      else if (strcmp(token,"endfacet") == 0) {
        assert(vertexVec.size()%3 == 0);
        const int iv0 = vertexVec.size()-3; assert(iv0 >= 0);
        // skip any collapsed tris...
        if ((vertexVec[iv0] == vertexVec[iv0+1])||(vertexVec[iv0] == vertexVec[iv0+2])||(vertexVec[iv0+1] == vertexVec[iv0+2])) {
          vertexVec.resize(iv0);
          assert(lastIvPtr);
          *lastIvPtr = iv0;
        }
        else {
          assert(int(edgeVec.size()) == iv0);
          edgeVec.resize(iv0+3);
          // look for edge pairs...
          for (int i = 0; i < 3; ++i) {
            map<const pair<DoubleVertex,DoubleVertex>,int>::iterator iter = edgeMap.find(pair<DoubleVertex,DoubleVertex>(vertexVec[iv0+i],vertexVec[iv0+(i+1)%3]));
            if (iter != edgeMap.end()) {
              assert(edgeVec[iter->second] == -1);
              edgeVec[iter->second] = iv0+i;
              edgeVec[iv0+i] = iter->second;
              // found a match...
              edgeMap.erase(iter);
            }
            else {
              // not found, so add the reverse edge to edgeMap...
              edgeMap[pair<DoubleVertex,DoubleVertex>(vertexVec[iv0+(i+1)%3],vertexVec[iv0+i])] = iv0+i;
              edgeVec[iv0+i] = -1;
            }
          }
        }
      }
      else if (strcmp(token,"vertex") == 0) {
        vertexVec.push_back(DoubleVertex());
        sscanf(getNextStlToken(fp,pos,max_pos,buf),"%lf",&vertexVec.back().x);
        sscanf(getNextStlToken(fp,pos,max_pos,buf),"%lf",&vertexVec.back().y);
        sscanf(getNextStlToken(fp,pos,max_pos,buf),"%lf",&vertexVec.back().z);
        if (!lastIvPtr) {
          // use stl filename for zone name...
          zoneNameMap[filename] = 0;
          newIzoneAndLastIvVec.push_back(pair<int,int>(0,vertexVec.size()));
          lastIvPtr = &(newIzoneAndLastIvVec.back().second);
        }
        *lastIvPtr = vertexVec.size();
      }
    }

    edgeMap.clear();
    fclose(fp);
    delete[] buf;

  }

  // check...

  assert(vertexVec.size()%3 == 0);
  const int nt = vertexVec.size()/3;
  assert(edgeVec.size() == vertexVec.size());

  cout << " > nt: " << nt << endl;

  // compress nodes...

  IntFlag flag(vertexVec.size());
  for (int iv = 0, nv_all = vertexVec.size(); iv < nv_all; ++iv) flag[iv] = iv;

  for (int ie = 0, ne = edgeVec.size(); ie < ne; ++ie) {
    if (edgeVec[ie] != -1) {
      const int ie_match = edgeVec[ie];
      assert(ie_match != ie);
      assert((ie_match >= 0)&&(ie_match < ne));
      assert(edgeVec[ie_match] == ie);
      if (ie_match > ie) {
        {
          int iv0 = ie;
          int iv1_match = ie_match+1; if (iv1_match%3 == 0) iv1_match -= 3;
          assert(vertexVec[iv0] == vertexVec[iv1_match]);
          int iv0_ = flag[iv0];
          while (iv0_ != flag[iv0_]) iv0_ = flag[iv0_];
          int iv1_match_ = flag[iv1_match];
          while (iv1_match_ != flag[iv1_match_]) iv1_match_ = flag[iv1_match_];
          flag[iv0_] = flag[iv1_match_] = min(iv0_,iv1_match_);
        }
        {
          int iv1 = ie+1; if (iv1%3 == 0) iv1 -= 3;
          int iv0_match = ie_match;
          assert(vertexVec[iv1] == vertexVec[iv0_match]);
          int iv1_ = flag[iv1];
          while (iv1_ != flag[iv1_]) iv1_ = flag[iv1_];
          int iv0_match_ = flag[iv0_match];
          while (iv0_match_ != flag[iv0_match_]) iv0_match_ = flag[iv0_match_];
          flag[iv1_] = flag[iv0_match_] = min(iv1_,iv0_match_);
        }
      }
    }
  }

  int nv = 0;
  for (int iv = 0, nv_all = vertexVec.size(); iv < nv_all; ++iv) {
    if (flag[iv] == iv) {
      ++nv;
      flag[iv] = -nv; // -1 indexing
    }
    else {
      int iv_ = flag[iv];
      while (iv_ >= 0) iv_ = flag[iv_];
      flag[iv] = iv_;
    }
  }

  // ========================================
  // set the simple surface stuff...
  // ========================================

  const int ss_nsp0 = nsp;

  nsp += nv;

  growNspData(nsp,ss_nsp0);

  nv = 0;
  for (int iv = 0, nv_all = vertexVec.size(); iv < nv_all; ++iv) {
    if (flag[iv] == -nv-1) {
      xsp[nv+ss_nsp0][0] = vertexVec[iv].x;
      xsp[nv+ss_nsp0][1] = vertexVec[iv].y;
      xsp[nv+ss_nsp0][2] = vertexVec[iv].z;
      ++nv;
    }
  }
  assert(nv+ss_nsp0 == nsp);
  vertexVec.clear();

  const int ss_nst0 = nst;
  nst += nt;

  growNstData(nst,ss_nst0);

  const int nzone = zoneNameMap.size();

  // cout << "nzone: " << nzone << endl;

  // start by adding all the zones to the ss...
  int * ss_zone_of_zone = new int[nzone];

  for (int jj = 0,jj_end=zoneVec.size(); jj < jj_end; ++jj) {
    map<const string,int>::iterator iter = zoneNameMap.find(zoneVec[jj].getName());
    if (zoneNameMap.find(zoneVec[jj].getName()) != zoneNameMap.end()) {
      // we found this zone already in the ss...
      ss_zone_of_zone[iter->second] = jj;
      zoneNameMap.erase(iter);
      //assert(0); // should be ok, but check // TODO decide what to do with duplicates
    }
  }
  // now everything left in the map is new. Add it in the sorted order...
  for (map<const string,int>::iterator iter = zoneNameMap.begin(); iter != zoneNameMap.end(); ++iter) {
    ss_zone_of_zone[iter->second] = zoneVec.size();
    zoneVec.push_back(SurfaceZone(iter->first));
  }

  int it_end = 0;
  for (int ii = 0,ii_end=newIzoneAndLastIvVec.size(); ii < ii_end; ++ii) {
    // recall that the terminal vertex index is stored in this Vec...
    const int it_begin = it_end;
    const int iv_end = newIzoneAndLastIvVec[ii].second;
    assert(iv_end%3 == 0);
    it_end = iv_end/3;
    assert(it_end > it_begin);
    const int izone = ss_zone_of_zone[newIzoneAndLastIvVec[ii].first];
    // now loop through the new tris associated with this zone...
    for (int it = it_begin; it != it_end; ++it) {
      FOR_I3 {
        spost[it+ss_nst0][i] = -flag[it*3+i]-1 + ss_nsp0;
      }
      znost[it+ss_nst0] = izone;
    }
  }

  DELETE(ss_zone_of_zone);
  return 0;
}

void SimpleSurface::writeSelectedTrisToStl(const string& filename, const bool single_zone) {
  if (mpi_rank == 0) {
    COUT1("writeSelectedTrisToStl(" << filename << ")");

    if (st_flag.count() == 0) {
      CWARN("no tris were flagged for export; skipping");
      return;
    }

    if (single_zone) COUT2(" > writing flagged tris to a single STL solid");
    char fname[128];
    sprintf(fname,"%s.stl",filename.c_str());
    FILE * fp = fopen(fname,"w");

    if (single_zone) writeSelectedTrisToStlSolid(fp, filename);
    else {
      // write flagged tris by zone to the file
      // first determine all the zones that need to be written
      zone_flag.resize(zoneVec.size());
      zone_flag.setAll(0);
      for (int ist=0; ist<nst; ++ist) {
        if (st_flag[ist]) {
          zone_flag[znost[ist]] = 1;
          st_flag[ist] = znost[ist];  // flag now holds which zone it is in
        }
        else st_flag[ist] = -1;  // don't use this tri
      }

      if (zone_flag.count()) {
        IntFlag zone_st_flag(nst);
        for (int izn=0,nzn=zoneVec.size(); izn<nzn; ++izn) {
          if (zone_flag[izn]) {
            zone_st_flag.setAll(0);
            for (int ist=0; ist<nst; ++ist) {
              if (st_flag[ist] == izn) zone_st_flag[ist] = 1;
            }
            assert(zone_st_flag.count());
            writeSelectedTrisToStlSolid(fp,zoneVec[izn].getName(),zone_st_flag);
          }
        }

      }
      else {
        CWARN("no tris were flagged for export; skipping");
      }
    }

    fclose(fp);
  }
}

void SimpleSurface::writeSelectedTrisToStlSolid(FILE * fp, const string& solidname, IntFlag& my_st_flag) const {
  if (mpi_rank == 0) {
    COUT1("writeSelectedTrisToStlSolid(" << solidname << ")");

    fprintf(fp,"solid %s\n",solidname.c_str());

    for (int ist = 0; ist < nst; ++ist) {
      if (my_st_flag[ist]) {
        fprintf(fp,"facet normal");
        double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
        FOR_I3 {
          const double mag_n = MAG(n);
          if (mag_n > 0.0) n[i] /= MAG(n);  // only scale normal if magnitude is valid
          fprintf(fp," %.6f",n[i]);
        }
        fprintf(fp,"\n");


        fprintf(fp,"  outer loop\n");
        FOR_I3 {
          fprintf(fp,"    vertex");
          int iv = spost[ist][i];
          FOR_J3 fprintf(fp," %.6f",xsp[iv][j]);
          fprintf(fp,"\n");
        }
        fprintf(fp,"  endloop\n");
        fprintf(fp,"endfacet\n");
      }
    }

    fprintf(fp,"endsolid %s\n",solidname.c_str());
  }
}

void SimpleSurface::writeSelectedTrisToStlSolid(FILE * fp, const string& solidname) const {
  if (mpi_rank == 0) {
    COUT1("writeSelectedTrisToStlSolid(" << solidname << ")");

    fprintf(fp,"solid %s\n",solidname.c_str());

    for (int ist = 0; ist < nst; ++ist) {
      if (st_flag[ist] != 0) {
        fprintf(fp,"facet normal");
        double n[3] = TRI_NORMAL_2(xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]);
        FOR_I3 {
          const double mag_n = MAG(n);
          if (mag_n > 0.0) n[i] /= MAG(n);  // only scale normal if magnitude is valid
          fprintf(fp," %.6f",n[i]);
        }
        fprintf(fp,"\n");


        fprintf(fp,"  outer loop\n");
        FOR_I3 {
          fprintf(fp,"    vertex");
          int iv = spost[ist][i];
          FOR_J3 fprintf(fp," %.6f",xsp[iv][j]);
          fprintf(fp,"\n");
        }
        fprintf(fp,"  endloop\n");
        fprintf(fp,"endfacet\n");
      }
    }

    fprintf(fp,"endsolid %s\n",solidname.c_str());
  }
}


#undef STL_HEADER_SIZE
#undef STL_SIZEOF_FACET
#undef STL_MIN_FILE_SIZE
#undef STL_LABEL_SIZE

#undef FAST_STL_BUFFER_SIZE
#undef FAST_STL_TOKEN_SIZE
