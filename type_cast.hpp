#pragma once

template<typename Parent>
auto&& parent_cast(auto&& child) {
    return static_cast<Parent>(
            const_cast<
                std::remove_cvref_t<decltype(child)>
            &>(
                child))
            ;
}
