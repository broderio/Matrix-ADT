//
//  MatrixErrors.h
//  Matrix
//
//  Created by Broderick Riopelle on 12/26/21.
//

#ifndef MatrixErrors_h
#define MatrixErrors_h

#include <stdio.h>
#include <string>

class dimension_mismatch {
private:
    std::string msg;
public:
    dimension_mismatch(std::string msg_in);
    std::string what() const;
};

class invalid_index {
private:
    std::string msg;
public:
    invalid_index(std::string msg_in);
    std::string what() const;
};

class rectangular_matrix {
private:
    std::string msg;
public:
    rectangular_matrix(std::string msg_in);
    std::string what() const;
};

class singular_matrix {
private:
    std::string msg;
public:
    singular_matrix(std::string msg_in);
    std::string what() const;
};

class empty_matrix {
private:
    std::string msg;
public:
    empty_matrix(std::string msg_in);
    std::string what() const;
};

dimension_mismatch::dimension_mismatch(std::string msg_in)
: msg(msg_in) {}

std::string dimension_mismatch::what() const
{
    return msg;
}

invalid_index::invalid_index(std::string msg_in)
: msg(msg_in) {}

std::string invalid_index::what() const
{
    return msg;
}

rectangular_matrix::rectangular_matrix(std::string msg_in)
: msg(msg_in) {}

std::string rectangular_matrix::what() const
{
    return msg;
}

singular_matrix::singular_matrix(std::string msg_in)
: msg(msg_in) {}

std::string singular_matrix::what() const
{
    return msg;
}

empty_matrix::empty_matrix(std::string msg_in)
: msg(msg_in) {}

std::string empty_matrix::what() const
{
    return msg;
}

#endif /* MatrixErrors_h */
