#pragma once

#include "Types.h"
#include <vector>

class ParticeList
{
public:
    ParticeList(void)
    {    };
    ~ParticeList(void)
    {    };

    struct SphereModel {
        openvdb::Vec3R pos;
        openvdb::Real  r;
    };

private:

    std::vector<SphereModel> m_ParticleList;

public:

    openvdb::Vec3R pos(int n)   const {return m_ParticleList[n].pos;}
    openvdb::Vec3R vel(int n)   const {return openvdb::Vec3R(1);}
    openvdb::Real radius(int n) const {return m_ParticleList[n].r;}
    
    void add(const openvdb::Vec3R pos, const openvdb::Real r)
    {
        SphereModel sp;
        sp.pos = pos;
        sp.r = r;
        m_ParticleList.push_back(sp);
    }
   

    size_t size() const { return m_ParticleList.size(); }
    void getPos(size_t n,  openvdb::Vec3R&pos) const { pos = m_ParticleList[n].pos; }

    
    void getPosRad(size_t n,  openvdb::Vec3R& pos, openvdb::Real& rad) const {
        pos = m_ParticleList[n].pos;
        rad = m_ParticleList[n].r;
    }
    void getPosRadVel(size_t n,  openvdb::Vec3R& pos, openvdb::Real& rad, openvdb::Vec3R& vel) const {
        pos = m_ParticleList[n].pos;
        rad = m_ParticleList[n].r;
        vel = openvdb::Vec3R(1);
    }
    // The method below is only required for attribute transfer
    void getAtt(size_t n, openvdb::Index32& att) const { att = n; }

};

