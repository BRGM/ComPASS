#pragma once

template <typename DestinationType, unsigned int parts = 1>
struct BufferInfo
{
    typedef DestinationType value_type;
    const value_type * const data;
    const std::size_t size;
    std::vector<std::size_t> shape;
    std::vector<std::size_t> stride;
    template <typename SourceType>
    BufferInfo(const std::vector<SourceType>& v) :
        data{ reinterpret_cast<const value_type *>(v.data()) },
        size{ sizeof(value_type) * parts * v.size() },
        shape{},
        stride{}
    {
        constexpr auto scalar_size = sizeof(value_type);
        static_assert(sizeof(SourceType) == parts * sizeof(value_type), "inconsistent memory size");
        static_assert(sizeof(SourceType) % scalar_size == 0, "inconsistent memory size");
        shape.push_back(v.size());
        stride.push_back(parts * scalar_size);
        if (parts > 1) {
            shape.push_back(parts);
            stride.push_back(scalar_size);
        }
    }
};
