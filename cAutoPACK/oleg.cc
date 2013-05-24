#include <vector>
#include <cstdlib>
#include <iostream>

struct big_grid { // needs 8*n bytes 
    std::vector<unsigned> all; 
    std::vector<unsigned> empty; 
    unsigned num_empty;
    big_grid(unsigned n) : all(n),
      empty(n),
      num_empty(n) {
    for(unsigned i = 0; i < n; ++i) { all[i] = i;
                empty[i] = i;
            }
    }

unsigned sample_empty() const {
    return empty[rand() % num_empty]; }

void set_filled(unsigned i) {
    unsigned empty_index = all[i];
    if(empty_index < num_empty) { // really is empty
        --num_empty;
        unsigned filled_index = empty[num_empty]; 
        std::swap(empty[empty_index], empty[num_empty]); 
        std::swap(all[i], all[filled_index]);
    } }

bool is_empty(unsigned i) const { 
    return all[i] < num_empty;
} };

template<typename T> void test(unsigned n) { T g(n);
for(unsigned i = 0; i < n; ++i) {
unsigned s = g.sample_empty(); g.set_filled(s);
if(i + 10 >= n) // print the last 10 indexes
            std::cout << i << ' ' << s << '\n';
} }
int main() {
const unsigned n = 100*100*100; test<small_grid>(n); // runs in 0.22 s test<big_grid>(n); // runs in 0.06 s
}
}
