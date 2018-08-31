#pragma once

template <typename Handle>
struct Functor_handle
{
    Handle handle;
    template <typename X>
    bool operator()(const X& x) const noexcept {
        return handle->operator(x);
    }
};