//
//  matrix_operators.cpp
//  Matrix
//
//  Created by Broderick Riopelle on 12/26/21.
//

#include <stdio.h>
#include <string>

class dimension_error {
private:
    std::string msg;
public:
    dimension_error(std::string msg_in);
    std::string what() const;
};

class index_error {
private:
    std::string msg;
public:
    index_error(std::string msg_in);
    std::string what() const;
};

class not_square_error {
private:
    std::string msg;
public:
    not_square_error(std::string msg_in);
    std::string what() const;
};

class is_singular_error {
private:
    std::string msg;
public:
    is_singular_error(std::string msg_in);
    std::string what() const;
};

class is_empty_error {
private:
    std::string msg;
public:
    is_empty_error(std::string msg_in);
    std::string what() const;
};

dimension_error::dimension_error(std::string msg_in)
: msg(msg_in) {}

std::string dimension_error::what() const
{
    return msg;
}

index_error::index_error(std::string msg_in)
: msg(msg_in) {}

std::string index_error::what() const
{
    return msg;
}

not_square_error::not_square_error(std::string msg_in)
: msg(msg_in) {}

std::string not_square_error::what() const
{
    return msg;
}

is_singular_error::is_singular_error(std::string msg_in)
: msg(msg_in) {}

std::string is_singular_error::what() const
{
    return msg;
}

is_empty_error::is_empty_error(std::string msg_in)
: msg(msg_in) {}

std::string is_empty_error::what() const
{
    return msg;
}

