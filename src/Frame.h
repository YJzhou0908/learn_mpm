#pragma once
#include "common.h"

class Frame {
public:
	int mWidth;
	int mHeight;
	int mChannels;

	Frame() {}

	Frame(int width, int height, int channels) {
		mWidth = width;
		mHeight = height;
		mChannels = channels;
		
	}

	void clear(){
		mData.clear();
	}

	void resize() {
		mData.resize(mWidth * mHeight * mChannels);
	}

	void drawColor(int x, int y, int c, float value) {
		if (x < 0 || x >= mWidth || y < 0 || y >= mHeight || c < 0 || c > 2)
			return;
		unsigned char color = (unsigned char)(value * 255);
		mData[y * mWidth * mChannels + x * mChannels + c] = color;

	}

	void outputImg(const std::string& filename) {
		stbi_flip_vertically_on_write(true);
		int result = stbi_write_jpg(filename.c_str(), mWidth, mHeight, mChannels, mData.data(), 90);
		return;
	}

	vector<unsigned char> mData;


};