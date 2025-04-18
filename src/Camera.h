#pragma once

#include "common.h"


using namespace glm;


class Camera {
public:
    dvec3 mPosition;
    dvec3 mLookAt;
    dvec3 mUp; 
    glm::dmat4 mProjectionViewMatrix;
    glm::dmat4 mViewMatrix;
    glm::dmat4 mProjectionMatrix;
    double mWidth;
    double mHeight;


public:
    Camera() {}

    Camera(glm::dvec3 position, glm::dvec3 lookAt, glm::dvec3 up) {
        mPosition = position;
        mViewMatrix = glm::lookAt(position, lookAt, up);
        mProjectionMatrix = glm::perspective<float>(glm::radians(45.0f), 1, 0.1f, 100.0f);

        mProjectionViewMatrix = mProjectionMatrix * mViewMatrix;

    }


    // return ndc position
    dvec3 project(const dvec3& point) {
       dvec4 ndc = mProjectionViewMatrix * glm::dvec4(point, 1.0) ;
       return dvec3(ndc / ndc.w);
    }
};
