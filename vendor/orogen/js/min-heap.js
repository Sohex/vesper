/**
 * Inline binary min-heap keyed on an external Float32Array of priorities.
 * Each cell is pushed/popped exactly once — no decrease-key needed.
 *
 * Shared by the drainage conditioning in terrain-post.js and the depression
 * detection in basins.js so both flood the terrain with identical ordering.
 * Deliberately has no tie-break on cell index: push order is a plain ascending
 * loop, so the result is already deterministic, and adding one would perturb
 * existing terrain output for no correctness gain.
 */
export class MinHeap {
    constructor(keyArray) {
        this._key = keyArray;
        this._data = [];
    }
    get size() { return this._data.length; }
    push(cell) {
        this._data.push(cell);
        let i = this._data.length - 1;
        while (i > 0) {
            const parent = (i - 1) >> 1;
            if (this._key[this._data[parent]] <= this._key[this._data[i]]) break;
            const tmp = this._data[parent]; this._data[parent] = this._data[i]; this._data[i] = tmp;
            i = parent;
        }
    }
    pop() {
        const top = this._data[0];
        const last = this._data.pop();
        if (this._data.length > 0) {
            this._data[0] = last;
            let i = 0;
            const n = this._data.length;
            for (;;) {
                const l = 2 * i + 1, r = 2 * i + 2;
                let smallest = i;
                if (l < n && this._key[this._data[l]] < this._key[this._data[smallest]]) smallest = l;
                if (r < n && this._key[this._data[r]] < this._key[this._data[smallest]]) smallest = r;
                if (smallest === i) break;
                const tmp = this._data[i]; this._data[i] = this._data[smallest]; this._data[smallest] = tmp;
                i = smallest;
            }
        }
        return top;
    }
}
