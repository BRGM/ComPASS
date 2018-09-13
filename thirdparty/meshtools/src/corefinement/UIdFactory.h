#pragma once

template <typename T>
struct UIdFactory
{
    typedef T Id_type;
    T count;
    UIdFactory() :
        count{ 0 }
    {}
    UIdFactory(const UIdFactory&) = delete;
    UIdFactory& operator=(const UIdFactory&) = delete;
    T make() {
        return count++;
    }
    constexpr T id_limit() const noexcept {
        return static_cast<T>(std::numeric_limits<T>::max());
    }
};

