//Functions to Operate on Vector Classes

#ifndef NL_VEC_FUNC_H
#define NL_VEC_FUNC_H

#include <iostream>
#include <cmath>

//Sum of all vector values
template <typename T>
double vSum(T vec){
    double sum=0.0;
    for(int i=0; i<vec.n; ++i) sum+=vec.v[i];
    return sum;
}

//Absolute sum of all vector values
template <typename T>
double vAsum(T vec){
    double asum=0.0;
    for(int i=0; i<vec.n; ++i) asum+=abs(vec.v[i]);
    return asum;
}

//Modulus of vector
template <typename T>
double vNorm(T vec){
    double norm;
    norm=sqrt(vSum(vec*vec));
    return norm;
}

//Square modulus of vector
template <class T>
double vNormSq(T vec){
    double normSq;
    normSq=vSum(vec*vec);
    return normSq;
}

//Mean of vector
template <typename T>
double vMean(T vec){
    return double(vSum(vec))/vec.n;
}

//Minimum value in vector
template <typename T>
double vMinimum(T vec){
    double min=vec[0];
    for(int i=1; i<vec.n; ++i){
        if(vec[i]<min) min=vec[i];
    }
    return min;
}

//Maximum value in vector
template <typename T>
double vMaximum(T vec){
    double max=vec[0];
    for(int i=1; i<vec.n; ++i){
        if(vec[i]>max) max=vec[i];
    }
    return max;
}

//Cyclic permutation of vector elements
template <typename T>
T vCyclicPermutation(T vec){
    T pVec(vec.n);
    for(int i=1; i<vec.n; ++i) pVec[i]=vec[i-1];
    pVec[0]=vec[vec.n-1];
    return pVec;
}

//Sort vector using heapsort algorithm (see numerical recipes)
template <typename T>
T vSort(T vec){
    int n = vec.n;
    if(n<=1) return vec;
    int i=n/2, parent, child;
    auto v=vec[0];
    T sVec(vec);
    if(n<2) return sVec;

    for(;;){
        if(i>0){
            i--;
            v=sVec[i];
        }
        else{
            n--;
            if(n==0) break;
            v=sVec[n];
            sVec[n]=sVec[0];
        }

        parent=i;
        child=i*2+1;
        while(child<n){
            if(child+1<n && sVec[child+1]>sVec[child]) ++child;
            if(sVec[child]>v){
                sVec[parent] = sVec[child];
                parent=child;
                child=parent*2+1;
            }
            else break;
        }
        sVec[parent]=v;
    }

    return sVec;
}

//Sort vector using heapsort algorithm (see numerical recipes) giving sort positions
template <typename T>
VecF<int> vArgSort(T vec, bool reverse=false){
    int n = vec.n;
    int i=n/2, parent, child;
    auto v=vec[0];
    int w=0;
    T sVec(vec);
    VecF<int> aVec(vec.n);
    for(int i=0; i<aVec.n; ++i) aVec[i]=i;
    if(n<2) return aVec;

    for(;;){
        if(i>0){
            i--;
            v=sVec[i];
            w=aVec[i];
        }
        else{
            n--;
            if(n==0) break;
            v=sVec[n];
            w=aVec[n];
            sVec[n]=sVec[0];
            aVec[n]=aVec[0];
        }

        parent=i;
        child=i*2+1;
        while(child<n){
            if(child+1<n && sVec[child+1]>sVec[child]) ++child;
            if(sVec[child]>v){
                sVec[parent] = sVec[child];
                aVec[parent] = aVec[child];
                parent=child;
                child=parent*2+1;
            }
            else break;
        }
        sVec[parent]=v;
        aVec[parent]=w;
    }
    if(reverse){
        for(int i=0, j=aVec.n-1; i<aVec.n; ++i,--j) sVec[i]=aVec[j];
        for(int i=0; i<aVec.n; ++i) aVec[i]=sVec[i];
    }
    return aVec;
}

//Find common values between vectors
template <typename T>
T vCommonValues(T vecA, T vecB){
    T vecC(vecA.n);
    int nCommon=0;
    for(int i=0; i<vecA.n; ++i){
        for(int j=0; j<vecB.n; ++j){
            if(vecA[i]==vecB[j]){
                vecC[nCommon]=vecA[i];
                ++nCommon;
            }
        }
    }
    T vecD(nCommon);
    for(int i=0; i<nCommon; ++i) vecD[i]=vecC[i];
    return vecD;
}

//Find if vector contains a value
template <typename T, typename S>
bool vContains(T vec, S val){
    bool contains=false;
    for(int i=0; i<vec.n; ++i){
        if(vec[i]==val){
            contains=true;
            break;
        }
    }
    return contains;
};

//Find unique values in vector
template <typename T>
T vUnique(T vecA){
    T vecB=vSort(vecA);
    int count=0;
    for(int i=1; i<vecB.n; ++i){
        if(vecB[i]==vecB[i-1]){
            vecB[i-1]=-1;
            ++count;
        }
    }
    T vecC(vecB.n-count);
    count=0;
    for(int i=0; i<vecB.n; ++i){
        if(vecB[i]>=0){
            vecC[count]=vecB[i];
            ++count;
        }
    }
    return vecC;
}

//Simple linear regression
template <typename T>
T vLinearRegression(T x, T y){
    if(x.n!=y.n) throw "Regression error - must be equal number of x and y values";
    T coefficients(3); //gradient, intercept and r-squared
    double sumX=0.0, sumY=0.0, sumXY=0.0, sumXX=0.0, sumYY=0.0;
    for(int i=0; i<x.n; ++i){
        sumX+=x[i];
        sumY+=y[i];
        sumXY+=x[i]*y[i];
        sumXX+=x[i]*x[i];
        sumYY+=y[i]*y[i];
    }
    double sqSumX=sumX*sumX;
    double sqSumY=sumY*sumY;
    double sdX=sqrt(x.n*sumXX-sqSumX);
    double sdY=sqrt(y.n*sumYY-sqSumY);
    double r=(x.n*sumXY-sumX*sumY)/(sdX*sdY);
    coefficients[0]=(r*sdY)/sdX;
    coefficients[1]=(sumY-coefficients[0]*sumX)/x.n;
    coefficients[2]=r*r;
    return coefficients;
}



#endif //NL_VEC_FUNC_H
