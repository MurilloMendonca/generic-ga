#pragma once

#include <string>
#include <cmath>
#include <algorithm>
namespace Helper{
	inline int gray2decimal(const std::string& gray){
		int decimal = 0;
		for(int i = 0; i < gray.size(); i++){
			decimal += (gray[i] - '0') * pow(2, gray.size() - i - 1);
		}
		return decimal;
	}

	inline std::string decimal2gray(int decimal, int size){
		std::string gray = "";
		for(int i = 0; i < size; i++){
			gray += (decimal % 2) + '0';
			decimal /= 2;
		}
		std::reverse(gray.begin(), gray.end());
		return gray;
	}

	inline std::string float2gray(float decimal, int size){
		std::string gray = "";
		for(int i = 0; i < size; i++){
			decimal *= 2;
			if(decimal >= 1){
				gray += '1';
				decimal -= 1;
			}
			else{
				gray += '0';
			}
		}
		return gray;
	}

	inline float gray2float(const std::string& gray){
		float decimal = 0;
		for(int i = 0; i < gray.size(); i++){
			decimal += (gray[i] - '0') * pow(2, -(i + 1));
		}
		return decimal;
	}
}