#pragma once
#include <cmath>
#include <algorithm>
#include <vector>
namespace Helper{

    inline int gray2decimal(const std::vector<bool>& gray){
        int decimal = 0;
        for(int i = 0; i < gray.size(); i++){
            decimal += gray[i] * pow(2, gray.size() - i - 1);
        }
        return decimal;
    }
    inline std::vector<bool> decimal2gray(int decimal, int numberOfBits){
        std::vector<bool> gray;
        for(int i = 0; i < numberOfBits; i++){
            gray.push_back((decimal >> i) ^ (decimal >> (i + 1)));
        }
        std::reverse(gray.begin(), gray.end());
        return gray;
    }


    inline std::vector<bool> float2gray(float decimal, int numberOfBits){
        std::vector<bool> gray;
        for(int i = 0; i < numberOfBits; i++){
            decimal *= 2;
            gray.push_back((int)decimal);
            decimal -= (int)decimal;
        }
        return gray;
    }

    inline float gray2float(const std::vector<bool>& gray){
        float decimal = 0;
        for(int i = 0; i < gray.size(); i++){
            decimal += gray[i] * pow(2, -(i + 1));
        }
        return decimal;
    }

}
