
#pragma once

#include <gestalt/sha2_256.h>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <iostream>

using namespace std;

#define SHA256_DIGEST_LENGTH 32

string getCurrentTimestamp() {
    time_t t = time(NULL);
    tm* timePtr = localtime(&t);

    char buf[20];
    snprintf(buf, sizeof(buf), "%04d%02d%02d%02d%02d%02d",
             1900 + timePtr->tm_year,
             1 + timePtr->tm_mon,
             timePtr->tm_mday,
             timePtr->tm_hour,
             timePtr->tm_min,
             timePtr->tm_sec);
    return string(buf);
}

// Hash file using Gestalt SHA256
int sha256_file(const char *path, char *outputBuffer) {
    if (!path || !outputBuffer) return -1;

    // Open the file
    ifstream file(path, ios::binary);
    if (!file.is_open()) {
        perror("Can't open file");
        return -1;
    }

    // Read entire file into a string
    stringstream buffer;
    buffer << file.rdbuf();
    string fileData = buffer.str();
    file.close();

    // Compute SHA256 using Gestalt
    string hash = hashSHA256(fileData);

    // Copy result to outputBuffer
    strncpy(outputBuffer, hash.c_str(), hash.size() + 1); // +1 for null terminator

    return 0;
}
