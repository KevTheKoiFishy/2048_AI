/** SHORTHAND REDEFINITIONS **/
const floor = Math.floor;
const round = Math.round;
const ceil  = Math.ceil;
const rand  = Math.random;
const sin   = Math.sin;
const π     = Math.PI;
const τ     = 2 * π;
const log   = console.log;
    const log_comb = log;
    const NNlog    = log;
    const gamerLog = log;
    const genLog   = log;
const table = console.table;

/** NATIVE PROTOTYPE METHOD DEFINITIONS **/
{
    /** NON-MUTATION */
    Number.prototype.abs = function() {
        return (this >= 0) ? this : -this;
    }

    Int8Array.prototype.equals = function(anotherArray) {
        if (this.length != anotherArray.length) return 0;
        for (let i = 0; i < this.length; ++i) {
            if (this[i] != anotherArray[i]) return 0;
        }
        return 1;
    }
    Int8Array.prototype.copy = function() {
        return this.map( (x) => x );
    }
    Int8Array.prototype.getGridString = function(n = 4) {
        let grid = this;

        let str = "";

        for (let y = 0; y < n; ++y) {
            for (let x = 0; x < n; ++x) {
                let power = grid[x + n*y];

                if (!power)  str += "     ";
                else         str += (1 << power).toString().padStart(5);
                if (x < n-1) str += " | ";
            }
            str += "\n";
        }

        return str;
    }

    Float32Array.prototype.copy = function() {
        return this.map( (x) => x );
    }
    Float32Array.prototype.equals = function(anotherArray) {
        if (this.length != anotherArray.length) return 0;
        for (let i = 0; i < this.length; ++i) {
            if (this[i] != anotherArray[i]) return 0;
        }
        return 1;
    }
    Float32Array.prototype.dot = function(anotherArray) {
        if (this.length != anotherArray.length) { throw "Dot product vector lengths don't match!"; }
        let sum = 0;
        this.forEach( (thisVal, index) => {
            if (thisVal == 0) return; // Allow prune via skipping
            sum += thisVal * anotherArray[index];
        });
        return sum;
    }
    Float32Array.prototype.multEach = function(anotherArray) {
        return this.map( (v, i) => v * anotherArray[i] );
    }
    Float32Array.prototype.plus = function(anotherArray) {
        return this.map( (v, i) => v + anotherArray[i] );
    }
    Float32Array.prototype.minus = function(anotherArray) {
        return this.map( (v, i) => v - anotherArray[i] );
    }
    Float32Array.prototype.scale = function(someScalar) {
        return this.map( (x) => x * someScalar );
    }
    Float32Array.prototype.SSE = function(targetOutput) {
        let sumOfSquares = 0;
        this.forEach( (x, i) => {
            const error = x - targetOutput[i];
            sumOfSquares += error * error;
        });
        return sumOfSquares;
    }

    Float32Array.prototype.max = function() {
        let maxVal = -Infinity;
        let maxInd = 0;
        this.forEach( (x, i) => {
            if (x < maxVal) return;
            maxVal = x;
            maxInd = i;
        });
        const out = [maxVal, maxInd];
        return out;
    }
    Float32Array.prototype.normSquare = function() {
        /**
         * Squares each array value and scales the values so the sum of squared elements is 1
         * This gives disproportionate weight to larger values
        */
        let magSquared = 0;

        let squaredArr = this.map( (x) => {
            if (x <= 0) return 0;
            let x_sq = x * x;
            magSquared += x_sq;
            return x_sq;
        });
        let normSquaredArr = squaredArr.scale((magSquared == 0) ? 0 : (1 / magSquared));
        return normSquaredArr;
    }
    Float32Array.prototype.randProbs = function() {
        /**
         * Takes normalized array as input and selects one of
         * the indexes by the probability distribution.
         * IE. [0.25, 0.25, 0.5].randProbs() returns
         * 0 with 25% chance
         * 1 with 25% chance
         * 2 with 50% chance
         */

        let percentile = rand();
        let probSum    = 0;
        let i          = 0;

        do { if (percentile < (probSum += this[i])) return i; }
        while (++i < this.length);

        return -1;
    }

    Float32Array.prototype.isZero = function(i = 0) {
        if      (this[i])              { return 0; }
        else if (i == this.length - 1) { return 1; }
        else                           { return this.isZero(i + 1); }
    }

    Float32Array.prototype.printStr = function(decimals = 4, scale = 1, unit = '') {
        return "[ " + Array.from(this, (v) => (v*scale).toFixed(decimals).toString() ).join(", ") + " ] " + unit;
    }

    /* MUTATION */
    Float32Array.prototype.plusEquals = function(anotherArray) {
        anotherArray.forEach( (x, i) => { this[i] += x });
    }
    Float32Array.prototype.minusEquals = function(anotherArray) {
        anotherArray.forEach( (x, i) => { this[i] -= x });
    }
    Float32Array.prototype.minusEquals_ignore0 = function(anotherArray) {   // Don't gradient-descent pruned weights
        anotherArray.forEach( (x, i) => { if(this[i] != 0) { this[i] -= x } });
    }
    Float32Array.prototype.scaleEquals = function(someScalar) {
        this.forEach( (x, i) => { this[i] *= someScalar } );
    }
    Float32Array.prototype.zeroOut = function() {
        this.forEach( (x, i) => { this[i] = 0 } );
    }

    Array.prototype.append = function(items) {
        if ( typeof(items) == "object" ) {
            items.forEach( (x) => { this.push(x) });
        } else {
            this.push(items);
        }
    }
    Array.prototype.appendAbs_ignore0 = function(items) {
        if ( typeof(items) == "object" ) {
            items.forEach( (x) => { if (x) { this.push(x.abs()) } });
        } else {
            if (items) { this.push(items.abs()) }
        }
    }
    Array.prototype.sliceStep = function(step = 1, start = 0, end = this.length) {
        let sliced = new Array();

        for (let i = start; i < end; i += step)
            sliced.push(this[i]);

        return sliced;
    }
}

/** GAME SIMULATION **/
class Game {
    constructor (n = 4, grid = new Int8Array(n * n).fill(0)) {
        this.n         = n;
        this.grid      = grid;
        this.score     = 0;
    }

    // DEBUG
    printGrid() {
        let n    = this.n;
        let grid = this.grid;

        let printGrid = Array.from( {length: n}, () => Uint32Array.from( {length: n}, () => 0) );

        for (let y = 0; y < n; ++y) {
            for (let x = 0; x < n; ++x) {
                let power = grid[x + n*y];
                if (!power) { printGrid[y][x] = " "; continue; }
                printGrid[y][x] = (1 << power);
            }
        }
        table(printGrid);
    }
    getGridString() {
        let n    = this.n;
        let grid = this.grid;

        let str = "";

        for (let y = 0; y < n; ++y) {
            for (let x = 0; x < n; ++x) {
                let power = grid[x + n*y];

                if (!power)  str += "     ";
                else         str += (1 << power).toString().padStart(5);
                if (x < n-1) str += " | ";
            }
            str += "\n";
        }

        return str;
    }

    // SIMULATION
    start() {
        this.score = 0;
        this.clearGrid();
        this.spawn();
    }
    clearGrid() {
        this.grid = new Int8Array(this.n * this.n).fill(0);
    }
    getFreeTiles() {
        this.freeTiles = new Array();
        this.grid.forEach( (V, I) => {
            if (!V) this.freeTiles.push(I);
        });
    }
    spawn() {
        let grid      = this.grid;

        let freeTiles = [];
        grid.forEach( (V, I) => {
            if (!V) freeTiles.push(I);
        });
        // *** log(this.getGridString());
        // *** log("freeTiles = ", freeTiles);
        let freeTileDiceRoll = floor(rand() * freeTiles.length);
        let freeTileToSpawn  = freeTiles[freeTileDiceRoll];
        
        if (rand() < 0.8) { grid[freeTileToSpawn] = 1; } // Spawn 2
        else              { grid[freeTileToSpawn] = 2; } // Spawn 4

        return freeTileToSpawn;
    }

    // Calls combine direction based on integers [0, 3] corresponding Up, Dwn, L, R
    combine(direction) {
        switch (direction) {
            case 0: this.combineUp();      break;
            case 1: this.combineDown();    break;
            case 2: this.combineLeft();    break;
            case 3: this.combineRight();   break;
        }
    }

    directionOf(direction) {
        switch (direction) {
            case 0: return "UP!";
            case 1: return "DOWN!";
            case 2: return "LEFT!";
            case 3: return "RIGHT!";
        }
    }

    combineUp() {
        let n = this.n;
        let grid = this.grid;

        // Handle each column
        for (let x = 0; x < n; ++x) {

            //Shift everything up
            // *** log_comb("  Start Shifting Up Column", x);
            let freeSpotAt = null;
            for (let y = 0; y < n; ++y) {
                let coord = y*n + x;        
                // *** log_comb("    Looking for free spot at", y, x, "coord =", coord);
                
                if (!grid[coord] && freeSpotAt == null) {  //Cannot use !freeSpotAt if freeSpotAt might be 0
                    freeSpotAt = coord;
                    // *** log_comb("    Found first free spot at", freeSpotAt);
                }
                else if ( grid[coord] && freeSpotAt != null) {
                    // *** log_comb("    Grid at", coord, "filled; free spot available.");
                    grid[freeSpotAt] = grid[coord];
                    // *** log_comb("    Setting free spot at", freeSpotAt, "to", grid[coord]);
                    grid[coord] = 0;
                    // *** log_comb("    Setting grid[", coord, "] to 0");
                    freeSpotAt  += n;
                    // *** log_comb("    Updating free spot position to 1 row below, at", freeSpotAt);
                }
                // *** this.printGrid();
            }
            // *** log_comb("  End Shifting Up Column", x);
            // *** this.printGrid();
            // *** log_comb();
            // *** log_comb();

            // Combine
            // *** log_comb("  Start Combining Column", x);
            for (let y = 0; y < n - 1; ++y) {
                let coord      = y*n + x;
                let coordBelow = coord + n;
                // *** log_comb("    Testing upwards combinability in cell (y,x) =", y, x);
                // *** log_comb("    coord =", coord, ", coordBelow =", coordBelow);

                // Skip combining if empty cell. Used to combine two 0s to make a 1 LMAO.
                if (!(grid[coord])) continue;
                // *** log_comb("    !(grid[coord]);", !(grid[coord]), "=> Not empty cell!");

                // Skip combining if different powers of 2
                if (grid[coord] != grid[coordBelow]) continue;
                // *** log_comb("    grid[coord] == grid[coordBelow];", grid[coord], "==", grid[coordBelow], "=> Combinable!");
                // *** this.printGrid();

                // Combine
                this.score += (1 << ++grid[coord]);
                // *** log_comb("    ++grid[coord]");
                // *** this.printGrid();

                // Starting with tile below augmented tile,
                // Replace each tile's value with value of tile below.
                // *** log_comb("    Moving tiles up!");
                for (let yMove = y + 1; yMove <= n - 2; ++yMove) {
                    let moveCoord      = yMove*n + x;
                    let moveCoordBelow = moveCoord + n;
                    // *** log_comb("      moveCoord =", moveCoord, ", moveCoordBelow =", moveCoordBelow);

                    // *** log_comb("      grid[moveCoord] = grid[moveCoordBelow];", grid[moveCoord], " <- ", grid[moveCoordBelow]);
                    grid[moveCoord] = grid[moveCoordBelow];
                    // *** this.printGrid();
                }
                // Bottom cell in this column becomes empty
                grid[(n-1)*n + x] = 0;
            }
            // *** log_comb("  Done Combining Column", x);
            // *** this.printGrid();
        }
    }
    combineDown() {
        let n = this.n;
        let grid = this.grid;

        // Handle each column
        for (let x = 0; x < n; ++x) {

            //Shift everything down
            let freeSpotAt = null;
            // *** log_comb("  Start Shifting Down Column", x);
            for (let y = n - 1; y >= 0; --y) {
                let coord = y*n + x;
                // *** log_comb("    Looking for free spot at", y, x, "coord =", coord);
                
                if (!grid[coord] && freeSpotAt == null) {
                    freeSpotAt = coord
                    // *** log_comb("    Found first free spot at", freeSpotAt);
                };
                if ( grid[coord] && freeSpotAt != null) {
                    // *** log_comb("    Grid at", coord, "filled; free spot available.");
                    grid[freeSpotAt] = grid[coord];
                    // *** log_comb("    Setting free spot at", freeSpotAt, "to", grid[coord]);
                    grid[coord] = 0;
                    // *** log_comb("    Setting grid[", coord, "] to 0");
                    freeSpotAt  -= n; // :)
                    // *** log_comb("    Updating free spot position to 1 row above, at", freeSpotAt);
                }
                // *** this.printGrid();
            }
            // *** log_comb("  End Shifting Down Column", x);
            // *** this.printGrid();
            // *** log_comb();
            // *** log_comb();

            // Combine
            // *** log_comb("  Start Combining Column", x);
            for (let y = n - 1; y >= 1; --y) {
                let coord      = y*n + x;
                let coordAbove = coord - n; // :)
                // *** log_comb("    Testing downwards combinability in cell (y,x) =", y, x);
                // *** log_comb("    coord =", coord, ", coordAbove =", coordAbove);

                // Skep empty cells
                if (!grid[coord]) continue; // :)
                // *** log_comb("    !(grid[coord]);", !(grid[coord]), "=> Not empty cell!");

                // Skip combining if different powers of 2
                if (grid[coord] != grid[coordAbove]) continue;
                // *** log_comb("    grid[coord] == grid[coordAbove];", grid[coord], "==", grid[coordAbove], "=> Combinable!");
                // *** this.printGrid();

                // Combine
                this.score += (1 << ++grid[coord]);
                // *** log_comb("    ++grid[coord]");
                // *** this.printGrid();

                // Starting with tile above augmented tile,
                // Replace each tile's value with value of tile above
                // *** log_comb("    Moving tiles Down!");
                for (let yMove = y - 1; yMove >= 1; --yMove) {
                    let moveCoord      = yMove*n + x;
                    let moveCoordAbove = moveCoord - n; // :)
                    // *** log_comb("      moveCoord =", moveCoord, ", moveCoordAbove =", moveCoordAbove);

                    // *** log_comb("      grid[moveCoord] = grid[moveCoordAbove];", grid[moveCoord], " <- ", grid[moveCoordAbove]);
                    grid[moveCoord] = grid[moveCoordAbove];
                    // *** this.printGrid();
                }
                // Top cell in this column becomes empty
                grid[x] = 0;

            }
        }
    }
    combineLeft() {
        let n = this.n;
        let grid = this.grid;

        // Handle each column
        for (let y = 0; y < n; ++y) {

            //Shift everything left
            // *** log_comb("  Start Shifting Left row", y);
            let freeSpotAt = null;
            for (let x = 0; x < n; ++x) {
                let coord = y*n + x;
                // *** log_comb("    Looking for free spot at", y, x, "coord =", coord);
                
                if (!grid[coord] && freeSpotAt == null) {  //Cannot use !freeSpotAt if freeSpotAt might be 0
                    freeSpotAt = coord;
                    // *** log_comb("    Found first free spot at", freeSpotAt);
                }
                else if ( grid[coord] && freeSpotAt != null) {
                    // *** log_comb("    Grid at", coord, "filled; free spot available.");
                    grid[freeSpotAt] = grid[coord];
                    // *** log_comb("    Setting free spot at", freeSpotAt, "to", grid[coord]);
                    grid[coord] = 0;
                    // *** log_comb("    Setting grid[", coord, "] to 0");
                    freeSpotAt  += 1;
                    // *** log_comb("    Updating free spot position to 1 cell right, at", freeSpotAt);
                }
                // *** this.printGrid();
            }
            // *** log_comb("  End Shifting left row", y);
            // *** this.printGrid();
            // *** log_comb();
            // *** log_comb();

            // Combine
            // *** log_comb("  Start Combining Row", y);
            for (let x = 0; x < n - 1; ++x) {        
                let coord      = y*n + x;
                let coordRight = coord + 1;
                // *** log_comb("    Testing leftwards combinability in cell (y,x) =", y, x);
                // *** log_comb("    coord =", coord, ", coordRight =", coordRight);

                // Skip combining if empty cell. Used to combine two 0s to make a 1 LMAO.
                if (!(grid[coord])) continue;
                // *** log_comb("    !(grid[coord]);", !(grid[coord]), "=> Not empty cell!");

                // Skip combining if different powers of 2
                if (grid[coord] != grid[coordRight]) continue;
                // *** log_comb("    grid[coord] == grid[coordRight];", grid[coord], "==", grid[coordRight], "=> Combinable!");
                // *** this.printGrid();

                // Combine
                this.score += (1 << ++grid[coord]);
                // *** log_comb("    ++grid[coord]");
                // *** this.printGrid();

                // Starting with tile below augmented tile,
                // Replace each tile's value with value of tile to right.
                // *** log_comb("    Moving tiles left!");
                for (let xMove = x + 1; xMove <= n - 2; ++xMove) {
                    let moveCoord      = y*n + xMove;
                    let moveCoordRight = moveCoord + 1;
                    // *** log_comb("      moveCoord =", moveCoord, ", moveCoordRight =", moveCoordRight);

                    // *** log_comb("      grid[moveCoord] = grid[moveCoordRight];", grid[moveCoord], " <- ", grid[moveCoordRight]);
                    grid[moveCoord] = grid[moveCoordRight];
                    // *** this.printGrid();
                }
                // Rightmost cell in this column becomes empty
                grid[y*n + (n-1)] = 0;
            }
            // *** log_comb("  Done Combining Row", y);
            // *** this.printGrid();
        }
    }
    combineRight() {
        let n = this.n;
        let grid = this.grid;

        // Handle each column
        for (let y = 0; y < n; ++y) {

            //Shift everything right
            // *** log_comb("  Start Shifting Right row", y);
            let freeSpotAt = null;
            for (let x = n - 1; x >= 0; --x) {
                let coord = y*n + x;
                // *** log_comb("    Looking for free spot at", y, x, "coord =", coord);
                
                if (!grid[coord] && freeSpotAt == null) {  //Cannot use !freeSpotAt if freeSpotAt might be 0
                    freeSpotAt = coord;
                    // *** log_comb("    Found first free spot at", freeSpotAt);
                }
                else if ( grid[coord] && freeSpotAt != null) {
                    // *** log_comb("    Grid at", coord, "filled; free spot available.");
                    grid[freeSpotAt] = grid[coord];
                    // *** log_comb("    Setting free spot at", freeSpotAt, "to", grid[coord]);
                    grid[coord] = 0;
                    // *** log_comb("    Setting grid[", coord, "] to 0");
                    freeSpotAt -= 1;
                    // *** log_comb("    Updating free spot position to 1 cell left, at", freeSpotAt);
                }
                // *** this.printGrid();
            }
            // *** log_comb("  End Shifting right row", y);
            // *** this.printGrid();
            // *** log_comb();
            // *** log_comb();

            // Combine
            // *** log_comb("  Start Combining Row", y);
            for (let x = n - 1; x >= 1; --x) {        
                let coord      = y*n + x;
                let coordLeft = coord - 1;
                // *** log_comb("    Testing leftwards combinability in cell (y,x) =", y, x);
                // *** log_comb("    coord =", coord, ", coordLeft =", coordLeft);

                // Skip combining if empty cell. Used to combine two 0s to make a 1 LMAO.
                if (!(grid[coord])) continue;
                // *** log_comb("    !(grid[coord]);", !(grid[coord]), "=> Not empty cell!");

                // Skip combining if different powers of 2
                if (grid[coord] != grid[coordLeft]) continue;
                // *** log_comb("    grid[coord] == grid[coordLeft];", grid[coord], "==", grid[coordLeft], "=> Combinable!");
                // *** this.printGrid();

                // Combine
                this.score += (1 << ++grid[coord]);
                // *** log_comb("    ++grid[coord]");
                // *** this.printGrid();

                // Starting with tile below augmented tile,
                // Replace each tile's value with value of tile to left.
                // *** log_comb("    Moving tiles right!");
                for (let xMove = x - 1; xMove >= 1; --xMove) {
                    let moveCoord      = y*n + xMove;
                    let moveCoordLeft  = moveCoord - 1;
                    // *** log_comb("      moveCoord =", moveCoord, ", moveCoordLeft =", moveCoordLeft);

                    // *** log_comb("      grid[moveCoord] = grid[moveCoordLeft];", grid[moveCoord], " <- ", grid[moveCoordLeft]);
                    grid[moveCoord] = grid[moveCoordLeft];
                    // *** this.printGrid();
                }
                // leftmost cell in this column becomes empty
                grid[y*n] = 0;
            }
            // *** log_comb("  Done Combining Row", y);
            // *** this.printGrid();
        }
    }

    // REMEMBER to update freeTiles[] upon combining
       // -> I feel like this would take more computation - updating freeTiles after each combine, shift, and spawn;
       // than counting all the free tiles after each move! Not doing it.
}

/** GFENETIC TRAINER **/
class Family {
    /**
     * If SEEDS provided, creates numCopies children from each seed
     * using NN.mitose() and NN.mutate() with successive span reduction.
     * Then picks the best children as the parents of the next generation.
     */

    constructor(numCopies = 10, numSurvivors = 5, seedGamers, numGames = 50) {
        this.seeds        = seedGamers;

        this.Igen         = 0;

        this.numCopies    = numCopies;
        this.numSurvivors = numSurvivors;        
        this.spanInit     = 0.1;
        this.spanRatio    = 0.1 ** ( 1 / (numCopies - 1) );

        this.numGames     = numGames;
        
        // Best NNPlayers
        this.Ψ = seedGamers || Array.from(
            {length: numSurvivors},
            (gamer, Igamer) => new NNGamer(0, Igamer)
        );
        // Children of best NNs
        this.φ = new Array(numSurvivors * numCopies);
    }

    stepGen() {
        let tick = performance.now();
        ++this.Igen;
        
        this.makeDaughters();
        this.playGames();
        this.selectFit();
        let tock = performance.now();
        genLog("Completed in", ((tock-tick)/1000).toFixed(3), "seconds.")
    }

    makeDaughters() {
        const Ψ = this.Ψ;
        const φ = this.φ;

        const numParents = Ψ.length;
        const numCopies  = this.numCopies;
        const spanInit   = this.spanInit;

        φ.length = 0; // Clear children

        for (let Imother = 0; Imother < numParents; ++Imother) {

            const thisMotherNN = Ψ[Imother].NN;
            let   span         = spanInit;

            for (let Idaughter = 0; Idaughter < numCopies; ++Idaughter) {

                const thisDaughter = thisMotherNN.mitose();
                thisDaughter.noiseProperty_uniform(span, 'α');
                thisDaughter.noiseProperty_uniform(span, 'bias');
                thisDaughter.noiseWeights_uniform(span);
                φ.push(new NNGamer(this.Igen, φ.length, undefined, thisDaughter));

                span *= this.spanRatio;

            }

        }

        // Include parents for comparison in future generations
        for (let Imother = 0; Imother < numParents; ++Imother) {
            φ.push(Ψ[Imother]);
        }

    }
    playGames() {
        const φ = this.φ;
        
        const Igen        = this.Igen;
        const numChildren = φ.length;
        const numGames    = this.numGames;

        genLog("GENERATION", Igen, "--------")
        const score_normalizer = 1 / this.numGames;

        for (let Igamer = 0; Igamer < numChildren; ++Igamer) {

            // *** genLog("  Igamer", Igamer, "Playing...")
            const gamerGal = φ[Igamer];

            // if (gamerGal.Igen < Igen && Igen > 1) continue; // Skip score-benchmarking the parents except the first time!

            // Track Best Game
            gamerGal.avgScore      = 0;
            gamerGal.avgEfficiency = 0;
            gamerGal.record        = 0;
            gamerGal.recordPtMoves = 0;
            gamerGal.bestGrid;

            for (let Igame = 0; Igame < numGames; ++Igame) {

                gamerGal.ready();                           // Clear grid + spawn tiles
                while (!gamerGal.makeMove(0)) {}            // Play until game over
                gamerGal.avgScore += gamerGal.game.score;   // Accumulate average score
                gamerGal.avgEfficiency += gamerGal.game.score / gamerGal.movesMade;
                
                if (gamerGal.game.score > gamerGal.record) { // Track Best Game - Just to see that the game is learning!
                    gamerGal.record        = gamerGal.game.score;
                    gamerGal.bestGrid      = gamerGal.game.grid.copy();
                    gamerGal.recordPtMoves = gamerGal.movesMade;
                }
                // *** genLog("    Game", Igame, "finished with score", gamerGal.game.score);
            }

            gamerGal.avgScore *= score_normalizer;
            gamerGal.avgEfficiency *= score_normalizer;
            // genLog("  Gamergal",s Igamer.toString().padStart(3), "averaged a score of", gamerGal.avgScore.toFixed(2).padStart(10), "with", gamerGal.avgEfficiency.toFixed(2).padStart(6), "pts/move over", numGames, "games");

        }

    }
    selectFit() {
        const Ψ = this.Ψ;
        const φ = this.φ;

        // φ.sort( (a, b) => b.avgScore - a.avgScore ); // Objective function = avg highscore
        φ.sort( (a, b) => b.avgEfficiency - a.avgEfficiency ); // Objective function = avg score / moves
        Ψ.length = 0;
        Ψ.append(φ.slice(0, this.numSurvivors));  // Best becomes parents for next generation

        genLog( "  GENERATION", this.Igen );
        genLog( "  HIGHSCORES: ", Array.from(φ, (gamerGal) => gamerGal.avgScore.toFixed(2) ) );
        genLog( "  BEST SCORE/MOVE: ", Array.from(φ, (gamerGal) => gamerGal.avgEfficiency.toFixed(2) ) );
        genLog( "  BEST GRID OF HIGHEST AVERAGING PLAYER: ", Ψ[0].record.toString().padStart(7), "pts with", Ψ[0].recordPtMoves.toString().padStart(5), "moves");
        genLog( Ψ[0].bestGrid.getGridString() );
    }
    
}

    /** ALGORITHMIC PLAYER - ON HOLD **/
        // Searches all possible game states and gets probabilistic
        // expected score for each key press after N key pressses
    class BFSGamer {
    }

    /** NEURAL PLAYER **/
    class NNGamer {
        constructor(Igen = 0, Ichild = 0, game_in = new Game(), NN_in = new NN()) {
            this.Igen          = Igen;
            this.Ichild        = Ichild;

            this.game          = game_in;
            this.NN            = NN_in;

            this.gameStates    = new Array();
            this.playerChoices = new Array();

            this.movesMade     = 0;
            
            this.avgScore      = 0;
            this.record        = 0;
            this.recordPtMoves = 0;
            this.bestGrid;
        }

        ready() {
            this.movesMade = 0;
            this.game.start();
        }

        makeMove (rememberMove = 1) {
            // *** gamerLog("");
            // *** gamerLog("--------------");
            // *** gamerLog("Current Score:", this.game.score);
            // *** gamerLog("Move number (0 indexed)", this.movesMade);
            
            let gameState0  = Float32Array.from(this.game.grid); // Type conversion to let it work with the NN
            let keyProbs    = this.NN.F(gameState0).copy();
            let keysTried   = 0b0000;
            let numKeysLeft = 4;
            let keyChoice;

            // *** gamerLog("Gamer sees:");
            // *** gamerLog(this.game.getGridString())
            // this.game.printGrid();
            // *** gamerLog("Neural Output:", keyProbs.printStr() )

            while (1) {
                keyProbs  = keyProbs.normSquare();

                keyChoice = keyProbs.randProbs();                 // Incentivezes networks to be confident in their answers.
                    // let [ , keyChoice]  = keyProbs.max();      // Causes player to get stuck if not combinable
                // *** gamerLog("Gamer's choices:", keyProbs.printStr(1, 100, "%") );
                // *** gamerLog("Gamer chooses", keyChoice, "=>", this.game.directionOf(keyChoice) );

                this.game.combine(keyChoice);
                keysTried |= (1 << keyChoice);
                
                // If something changed, then we made a successful move. Next move!
                if (!this.game.grid.equals(gameState0)) {
                    // *** gamerLog("  Move was successful! Onto next move.");
                    this.game.spawn();
                    break;
                }
                // If nothing changed, try other keys! If none causes a change, it's game over.
                else {
                    // If the current key choice made no moves, get rid of it and try again
                    keyProbs[keyChoice] = 0;
                    --numKeysLeft;

                    // *** gamerLog("  Move unsuccessful! keysTried:", keysTried.toString(2).padStart(4, "0"), "| numKeysLeft:", numKeysLeft);
                    
                    // If all keys have been tried with no luck, GAME OVER!
                    if (!numKeysLeft) {
                        // *** gamerLog("  No keys left, GAME OVER!");
                        this.handleDeath();
                        return 1;
                    }
                    // If not all keys tried, but probs are all zero, assign each remaining key same probability
                    if (keyProbs.isZero()) {
                        let keyProb = 1 / numKeysLeft;
                        // *** gamerLog("  Not all keys tried, but probs are all zero! Assign each remaining key to", (keyProb*100).toFixed(1), "%");
                        keyProbs.forEach ( (v, i) => {
                            if (keysTried & (1 << i)) return;
                            keyProbs[i] = keyProb;
                        });
                    }
                    // If not all keys tried, and there's >= 1 positive probability, go again
                }
            }

            ++this.movesMade;

            if (!rememberMove) return 0;
            this.gameStates.push(gameState0);
            this.playerChoices.push(keyChoice);

            return 0;
        }

        handleDeath () {
            // *** gamerLog("Ooops we died!");
            // If multithreading, onto next.

            return 0;
        }
        
    }

        /** NETWORK PRUNER **/
        class Gardener {
            constructor(NNToPrune = new NN(), inputs, outputs) {
                this.NN         = NNToPrune;
                this.inputs     = inputs;
                this.outputs    = outputs;
            }

            // Percentile = Percentage to remove.
            pruneWeights_magnitude(percentile) {
                /**
                 *  Gets absolute value of weights into a sorted array,
                    then finds a threshold weight based on percentile
                    such that weights whose are < thres will be eliminated
                    (set to 0).
                */

                // Gather all weights
                const NN        = this.NN;
                const net       = NN.net;

                const weights   = new Array();

                const numLayers = NN.numLayers;

                for (let Ilayer = 1; Ilayer < numLayers; ++Ilayer) {
                    
                    const thisLayer = net[Ilayer];
                    const numNodesThisLayer = NN.nodesByLayer[Ilayer];
                    const numNodesPrevLayer = NN.nodesByLayer[Ilayer - 1];
    
                    for (let Inode = 0; Inode < numNodesThisLayer; ++Inode) {
    
                        const thisNode = thisLayer[Inode];
                        weights.appendAbs_ignore0( thisNode.weights.scale(numNodesPrevLayer) ); // Undo initialization normalization
                        
                    }
                }

                // Sort weights
                weights.sort((a, b) => a - b);
                
                // Pick threshold weight
                const weightThres = weights[floor(weights.length * percentile)];

                // Sever connections below threshold
                for (let Ilayer = 1; Ilayer < NN.numLayers; ++Ilayer) {
                    
                    const thisLayer = net[Ilayer];
                    const numNodesThisLayer = NN.nodesByLayer[Ilayer];
                    const numNodesPrevLayer = NN.nodesByLayer[Ilayer - 1];

                    for (let Inode = 0; Inode < numNodesThisLayer; ++Inode) {

                        const thisNode = thisLayer[Inode];
                        const theseWeights = thisNode.weights;  // remember to unnorm here too
                        
                        for (let Iweight = 0; Iweight < numNodesPrevLayer; ++Iweight) {
                            const thisWeight = theseWeights[Iweight];
                            // *** NNlog("thisWeight =", thisWeight, " | weightThres =", weightThres);
                            
                            // Skip already pruned weights
                            if (thisWeight == 0)
                                continue;
                            // Delete weights who didn't make cutoff
                            if ( (thisWeight*numNodesPrevLayer).abs() < weightThres )
                                theseWeights[Iweight] = 0;
                        }
                    }
                }

            }
            pruneNodes_Ᾱ(percentile) {
                /**
                 * Gets the average absolute output of each neuron
                 * across all training data into a sorted array.
                 * Keeps neurons with top percentile of activation.
                 */
                NNlog("  Pruning Nodes by Ᾱ");

                const NN      = this.NN;
                const net     = NN.net;
                const inputs  = this.inputs;

                // MEASURE ACTIVATIONS OF EACH NEURON.
                const numSamples = inputs.length;
                const sample_normalizer = 1 / numSamples;

                // Init average activation to 0
                NN.clearProperty("Ᾱ");

                // Calculate activations
                for (let Isample = 0; Isample < numSamples; ++Isample) {

                    const thisInput = inputs[Isample];
                    NN.F(thisInput);``
                    
                    for (let Ilayer = 1; Ilayer < NN.numLayers - 1; ++Ilayer) {
                    
                        const thisLayer = net[Ilayer];
                        const numNodesThisLayer = NN.nodesByLayer[Ilayer];

                        for (let Inode = 0; Inode < numNodesThisLayer; ++Inode) {

                            const thisNode = thisLayer[Inode];
                            thisNode.Ᾱ += thisNode.out.abs();
                            
                        }
                    }
                }
                NN.scaleProperty("Ᾱ", sample_normalizer);

                // PRUNING
                // Rank avg activations
                const Ᾱs = NN.gatherPropertyAbs_ignore0("Ᾱ").sort((a, b) => a - b);
                NNlog("  Average Activations Before Deleting: ", Ᾱs);
                
                // Get cutoff
                const Ᾱ_thres = Ᾱs[round(Ᾱs.length * percentile)];
                NNlog("  Average Activations Threshold (", (percentile*100).toFixed(1), "% ): ", Ᾱ_thres);

                // Delete nodes less than cutoff
                this.pruneNodes_abs("Ᾱ", Ᾱ_thres);
                NNlog("  Average Activations After Deleting: ", NN.gatherPropertyAbs_ignore0("Ᾱ"));

            }
            pruneNodes_ΔC(percentile) {
                /**
                 * Gets the average change in cost of removing
                 each neuron, across all training data.
                * Keeps those with the highest cost impact.
                */
                NNlog("  Pruning Nodes by ΔC");

                const NN      = this.NN;
                const net     = NN.net;
                const inputs  = this.inputs;

                // MEASURE COST IMPACT OF EACH NEURON.
                const numSamples = inputs.length;
                const sample_normalizer = 1 / numSamples;

                // Init Cost Impact to 0
                NN.clearProperty("ΔC");

                // Calculate cost impact data
                for (let Isample = 0; Isample < numSamples; ++Isample) {

                    const thisInput = inputs[Isample];
                    const originalOutput = NN.F(thisInput).copy();
                    
                    for (let Ilayer = 1; Ilayer < NN.numLayers - 1; ++Ilayer) {
                    
                        const thisLayer = net[Ilayer];
                        const numNodesThisLayer = NN.nodesByLayer[Ilayer];

                        for (let Inode = 0; Inode < numNodesThisLayer; ++Inode) {

                            const thisNode = thisLayer[Inode];

                            const originalWeights = thisNode.weights.copy();
                            const originalBias    = thisNode.bias;
                            
                            // Pretend this neuron doesn't exist
                            thisNode.weights.zeroOut();
                            thisNode.bias = 0;
                            
                            // Get network output
                            const newOutput = NN.F(thisInput);
                            
                            // Calculate change in error
                            /*   
                                | originalOutput.SSE(originalOutput) - netOut.SSE(OriginalOutput) |
                                = | 0 - netOut.SSE(OriginalOutput) |
                                = | netOut.SSE(OriginalOutput) |
                            */
                            const inducedError   = newOutput.SSE(originalOutput);
                            thisNode.ΔC += inducedError;    // Accumulate the cost impact of removing each neuron
                            
                            // Put back original parameters
                            thisNode.weights = originalWeights;
                            thisNode.bias    = originalBias;
                        }
                    }
                }
                NN.scaleProperty("ΔC", sample_normalizer);

                // PRUNING
                // Rank cost impacts
                const ΔCs = NN.gatherPropertyAbs_ignore0("ΔC").sort((a, b) => a - b);
                NNlog("  Average Cost Impact Before Deleting: ", ΔCs);

                // Get cutoff
                const ΔC_thres = ΔCs[floor(ΔCs.length * percentile)];
                NNlog("  Average Cost Impact Threshold (", (percentile*100).toFixed(1), "% ): ", ΔC_thres);

                // Delete nodes less than cutoff
                this.pruneNodes_abs("ΔC", ΔC_thres);
                NNlog("  Average Cost Impact After Deleting: ", NN.gatherPropertyAbs_ignore0("ΔC"));
            }

            pruneNodes_abs(parameter = "ΔC", threshold = 0.1) {
                const NN      = this.NN;
                const net     = NN.net;

                // Delete nodes less than cutoff
                    // Fixed error: Prevent pruning of output layer
                for (let Ilayer = 1; Ilayer < NN.numLayers - 1; ++Ilayer) {
                    
                    const thisLayer = net[Ilayer];

                    for (let Inode = 0; Inode < NN.nodesByLayer[Ilayer]; ++Inode) {

                        const thisNode = thisLayer[Inode];
                        if (thisNode[parameter].abs() < threshold) {
                            NN.removeNode(thisNode.Ilayer, thisNode.Inode);
                            --Inode;
                        }
                        
                    }
                }
            }
            
        }

        /** SUPERVISED SGD HANDLER **/
        class Teacher {
            constructor
            (
                NNToTrain         = new NN(),
                inputs            = [new Float32Array()],
                outputs           = [new Float32Array()],
                outputTokens      = new Array(),
                learnRate         = 0.01,
                batchSize         = 5,
                numThreads        = 8
            )
            {
                this.NN           = NNToTrain;
                this.inputs       = inputs;
                this.outputs      = outputs;
                this.outputTokens = outputTokens;
                
                this.numSamples   = inputs.length;

                this.sIn          = new Array(this.numSamples);
                this.sOut         = new Array(this.numSamples);

                this.NN.learnRate = learnRate;
            }
            
            setVectorOutputs(mapFxn = function(outputToken) { return new Float32Array(10); }) {
                /**
                  * Turns tokens into corresponding vectors and adds
                    them to this.outputs using a user specified fxn.

                  * Sample mapFxn:
                    (outputToken) => {
                        let outputVect = new Float32Array(10);
                            outputVect[outputToken] = 1;
                        return outputVect;
                    }

                    ex. token 1 becomes [0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
                        token 2 becomes [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
                        etc.
                 */

                const outputs       = this.outputs;
                const outputTokens  = this.outputTokens;

                outputTokens.forEach( (thisToken, Isample) => {
                    outputs[Isample] = mapFxn(thisToken);
                });
            }

            SGD(numEpochs = 10, targetCost = 0, targetAccuracy = 1, toLog = true, tokenizeFxn = (outputVector) => outputVector.max()[1] ) {
                
                const NN    = this.NN;
                const sIn   = this.sIn;
                const sOut  = this.sOut;

                let log_local = log;
                if (!toLog) log_local = function() {};
                
                let tick = 0, tock = 0;
                let Iepoch = 0;

                for (Iepoch = 0; Iepoch < numEpochs; ++Iepoch) {

                    const [cost, accuracy] = this.getCostAccuracy( tokenizeFxn );
                    NN.cost = cost;
                    NN.accuracy = accuracy;
                    NNlog("Epochs Elapsed:", Iepoch, "(", (tock - tick).toFixed(3), "ms )", "Cost", cost.toFixed(5), "Accuracy", (accuracy * 100).toFixed(4), "%");
                    if (cost     < targetCost)     break;
                    if (accuracy > targetAccuracy) break;

                    this.scrambleTraining();
                    
                    tick = performance.now();

                    sIn.forEach( (thisInputVect, Idatum) => {
                        const wantedOutputVect = sOut[Idatum];

                        NN.F(thisInputVect);
                        NN.G(wantedOutputVect);
                        NN.D();
                        // NNlog(`At Iepoch ${Iepoch.toString().padStart(3)}, Idatum ${Idatum.toString().padStart(4)}, cost is ${myNN.cost.toPrecision(5).padStart(12)}, output =`,
                        //     Array.from(myNN.netOut, (x) => x.toPrecision(3).toString().padEnd(12) ).join("")
                        // );
                    } );

                    tock = performance.now();
                }

                const [cost, accuracy] = this.getCostAccuracy( tokenizeFxn );
                NN.cost = cost;
                NN.accuracy = accuracy;
                log_local("Epochs Elapsed:", Iepoch, "(", (tock - tick).toFixed(3), "ms )", "Cost", cost.toFixed(5), "Accuracy", (accuracy * 100).toFixed(4), "%");
            
            }
            scrambleTraining() {
                /**
                 * References this.sIn and this.sOut's elements 
                   to the original array's elements in random order.
                 */

                const numSamples    = this.numSamples;

                // Get scrambled order
                const unusedIndexes = Array.from({length: numSamples}, (x, i) => i);
                const trainingOrder = new Uint16Array(numSamples);

                for (let Idatum = 0; Idatum < numSamples; ++Idatum) {
                    const pickedIndexIndex      = floor(unusedIndexes.length * rand());
                    const pickedTrainingIndex   = unusedIndexes[ pickedIndexIndex ];
                    trainingOrder[Idatum]       = pickedTrainingIndex;
                    unusedIndexes.splice(pickedIndexIndex, 1);
                }

                // Make scrambled input and output array that
                // references the original arrays in specified order
                const In   = this.inputs;
                const Out  = this.outputs;
                const sIn  = this.sIn;
                const sOut = this.sOut;

                for (let Idatum = 0; Idatum < numSamples; ++Idatum) {
                    sIn[Idatum]  = In[ trainingOrder[Idatum] ];
                    sOut[Idatum] = Out[ trainingOrder[Idatum] ];
                }
            }
            getCostAccuracy(tokenizeFxn = function(outputVector = new Float32Array(10)) { return 0; }) {
                /**
                  * Gets performance of neural net in terms of
                    cost and accuracy.
                  * Cost = sample-averaged SSE between desired
                    and measured output layer
                  * Accuracy = compares output token, which is
                    determined by some objective function
                    tokenizeFxn(outputVector), with the given
                    outputTokens.
                 */
                const NN                 = this.NN;

                const numSamples         = this.numSamples;
                const inputs             = this.inputs;
                const wantedOutputTokens = this.outputTokens;

                let cost     = 0;
                let accuracy = 0;

                wantedOutputTokens.forEach( (wantedOutputToken, Isample) => {
                    let wantedOutputVect = this.outputs[Isample];

                    let realOutputVect  = NN.F(inputs[Isample]);
                    let realOutputToken = tokenizeFxn(realOutputVect) || realOutputVect.max()[1];

                    cost     +=  wantedOutputVect.SSE(realOutputVect);
                    accuracy += (wantedOutputToken == realOutputToken);
                });

                cost     /= numSamples;
                accuracy /= numSamples;
                
                return [cost, accuracy];
            }

        }

        /** NEURAL NET WITH VERSATILE SET OF METHODS **/
        class NN {
            constructor(nodesByLayer = [16, 16, 8, 4], toInitialize = true) {
                this.nodesByLayer   = nodesByLayer;
                this.numLayers      = nodesByLayer.length;
                this.outputsByLayer = Array.from( nodesByLayer, (numNodesThisLayer) => {new Float32Array(numNodesThisLayer)} );
                this.netOut         = new Float32Array(nodesByLayer[this.numLayers - 1]);

                this.cost           = 0;
                this.accuracy       = 0;
                this.learnRate      = 0.01;

                // Ever since I discovered these Array methods I've all but forgone for loops 😆
                this.net = Array.from(
                    nodesByLayer,
                    (numNodesThisLayer, Ilayer) => {
                        const thisLayer =
                            (Ilayer == 0)
                            ? new Float32Array(numNodesThisLayer)
                            : new Array(numNodesThisLayer).fill(0).map(
                                (thisNode, Inode) => new Node(Ilayer, Inode, nodesByLayer[Ilayer - 1], undefined, undefined, toInitialize)
                            );
                        return thisLayer;
                    }
                );
            }
            
            // BASIC FUNCTIONS
                forEachNode(fxn = (thisNode = new Node(), Ilayer = 1, Inode = 0) => {}) {
                    const net = this.net;
                    for (let Ilayer = 1; Ilayer < net.length; ++Ilayer) {
                        const thisLayer = net[Ilayer];

                        for (let Inode = 0; Inode < thisLayer.length; ++Inode) {
                            const thisNode = thisLayer[Inode];
                            fxn(thisNode, Ilayer, Inode);
                        }
                    }
                }

            // GD FUNCTIONS
                // Forward prop network
                F(netInput = new Float32Array(this.nodesByLayer[0])) {
                    // *** NNlog(netInput);
                    this.net.forEach( (thisLayer, Ilayer) => {
                        thisLayer = this.net[Ilayer];

                        // Store input layer
                        if (Ilayer == 0) {
                            this.net[0] = netInput.copy();
                            this.outputsByLayer[0] = this.net[0];
                        }
                        // Run through hidden layers
                        else {
                            const layerInput = (Ilayer > 1)
                                ? this.outputsByLayer[Ilayer - 1]
                                : this.net[0];

                            thisLayer.forEach( (thisNode, Inode) => {
                                thisNode = thisLayer[Inode];
                                thisNode.F(layerInput);
                            });

                            this.outputsByLayer[Ilayer] = this.getLayerOutput(Ilayer);
                            // *** NNlog("Layer", Ilayer, "has output =", this.outputsByLayer[Ilayer]);
                        }
                    });
                    
                    // Put output in array
                    this.netOut = this.outputsByLayer[this.numLayers - 1];
                    return this.netOut;
                }
                // Get gradient. Call F first!
                G(outputWeWant = new Float32Array(this.netOut.length).fill(1)) {
                    this.update_cost(outputWeWant);
                    this.update_output_dC_dA(outputWeWant);
                    this.update_G_hidden();
                }
                // Gradient Descent. Call F and G first!
                D() {
                    for (let Ilayer = this.numLayers - 1; Ilayer > 0; --Ilayer) {
                        const thisLayer = this.net[Ilayer];
                        
                        thisLayer.forEach( (thisNode, Inode) => {
                            thisNode = thisLayer[Inode];
                            thisNode.D();
                        });
                    }
                }
            
            // FORWARD PROP FUNCTIONS
                // convert node().out s to an array
                getLayerOutput(Ilayer) {
                    let outs = Float32Array.from(
                        this.net[Ilayer],
                        (thisNode) => thisNode.out
                    );
                    return outs;
                }

            // BACKPROP FUNCTIONS
                // Calculate cost and update last layer's dC_dA
                update_cost(outputWeWant = new Float32Array()) {
                    this.cost = this.netOut.SSE(outputWeWant);
                }
                update_output_dC_dA(outputWeWant = new Float32Array()) {
                    const outLayer = this.net[this.numLayers - 1];

                    outLayer.forEach( (thisNode, Inode) => {
                        thisNode = outLayer[Inode];
                        thisNode.dC_dA = this.learnRate * (thisNode.out - outputWeWant[Inode]);
                        /**
                         * d/dA( (A - A*)^2 ) = 2 (A - A*), where A is current output and A* is desired output.
                            * if A > A*, dC_dA is +, so subtracting gradient will bring A LOWER  to A*
                            * if A < A*, dC_dA is -, so subtracting gradient will bring A HIGHER to A*
                            * I omit the 2 to save computation. Just change learning rate LOL.
                        
                        * Okay I just realized that it'd be a lot more efficient to multiply by
                        learn rate here and letting it propagate through the network,
                        since all gradients will be proportional to dC_dA of final layer.
                            * Unless learnRate is 1, the gradient values will NOT be the real gradient,
                            rather the Δ gradient to add during gradient descent...
                            AHHHHHH PLEASE DON'T KILL ME FOR BAD NOTATION
                        */
                    });
                }
                // Each layer gets its gradient updated except the preceding layer's dC_dA except first layer
                update_G_hidden() {
                    for (let Ilayer = this.numLayers - 1; Ilayer > 0; --Ilayer) {
                        const thisLayer         = this.net[Ilayer];
                        const thisLayerIn       = this.outputsByLayer[Ilayer - 1];
                        const numNodesThisLayer = thisLayer.length;
                        
                        for  (let Inode = 0; Inode < numNodesThisLayer; ++Inode) {
                            const thisNode   = thisLayer[Inode];
                            const this_dC_dI = thisNode.G(thisLayerIn);
                            // ... Except first hidden layer!
                            if (Ilayer > 1) this.update_prevLay_dC_dAs(Ilayer, this_dC_dI);
                        }
                    }
                }
                // To the previous layer, add up how much the current layer's nodes "want" the previous layer's outputs to change
                update_prevLay_dC_dAs(Ilayer, stuffToAdd = new Float32Array()) {
                    const thisLayer = this.net[Ilayer];
                    const prevLayer = this.net[Ilayer - 1];

                    // Init to 0
                    prevLayer.forEach( (thisNode, InodePrevLay) => {
                        thisNode = prevLayer[InodePrevLay];
                        thisNode.dC_dA = 0;
                    });

                    // Accumulate
                    thisLayer.forEach( (thisNode, InodeThisLay) => {
                        thisNode = thisLayer[InodeThisLay];
                        const changesToPrev = thisNode.dC_dI;
                        changesToPrev.forEach( (delta_dC_dA, InodePrevLay) => {
                            try {
                                const nodeToChange = prevLayer[InodePrevLay];
                                nodeToChange.dC_dA += delta_dC_dA;
                            }
                            catch (err) {
                                NNlog(Ilayer, InodeThisLay);
                            }
                        });
                    });
                }
                
            // PRUNING & STAT FUNCTIONS
                removeNode(Ilayer = 1, Inode_rm = 0) {
                    // Remove this neuron from this layer
                    const net = this.net;

                    const changedLayer = net[Ilayer];
                          changedLayer.splice(Inode_rm, 1);

                    // Update nodesByLayer
                    --this.nodesByLayer[Ilayer];

                    // Update successive neuron's Inode property
                    let numNodesChangedLayer = changedLayer.length; // Will be one less than when the fxn was called
                    for (let Inode_updateID = Inode_rm; Inode_updateID < numNodesChangedLayer; ++Inode_updateID) {
                        const thisNode = changedLayer[Inode_updateID];
                        --thisNode.Inode;
                    }

                    // Remove this neuron's weights and gradient elements from next layer
                        // Skip if last layer
                        if (Ilayer == this.numLayers - 1) return;

                    const nextLayer = net[Ilayer + 1];
                    nextLayer.forEach( (thisNode, Inode_rmWeights) => {
                        thisNode = nextLayer[Inode_rmWeights];
                        
                        const new_Weights = Array.from(thisNode.weights);
                              new_Weights.splice(Inode_rm, 1)
                        thisNode.weights = Float32Array.from(new_Weights);

                        const new_dC_dI = Array.from(thisNode.dC_dI);
                              new_dC_dI.splice(Inode_rm, 1)
                        thisNode.dC_dI = Float32Array.from(new_dC_dI);

                        const new_dC_dW = Array.from(thisNode.dC_dW);
                              new_dC_dW.splice(Inode_rm, 1)
                        thisNode.dC_dW = Float32Array.from(new_dC_dW);

                        const new_dZ_dW = Array.from(thisNode.dZ_dW);
                              new_dZ_dW.splice(Inode_rm, 1)
                        thisNode.dZ_dW = Float32Array.from(new_dZ_dW);
                    });
                }
                clearProperty(propertyName = "") {
                    /**
                     * Sets a given numerical property's value to 0 for all nodes.
                     * Usually this property is used as an accumulator.
                    **/
                    const net = this.net;
                    const nodesByLayer = this.nodesByLayer;

                    for (let Ilayer = 1; Ilayer < net.numLayers; ++Ilayer) {
                    
                        const thisLayer = net[Ilayer];
                        const numNodesThisLayer = nodesByLayer[Ilayer];
                        
                        for (let Inode = 0; Inode < numNodesThisLayer; ++Inode) {
        
                            const thisNode = thisLayer[Inode];
                            thisNode[propertyName] = 0;
                            
                        }
                    }
                }
                scaleProperty(propertyName = "", scaleFactor = 0.01) {
                    const numLayers     = this.numLayers;
                    const nodesByLayer  = this.nodesByLayer;
                    const net           = this.net;

                    for (let Ilayer = 1; Ilayer < numLayers; ++Ilayer) {
                        const thisLayer = net[Ilayer];
                        const numNodesThisLayer = nodesByLayer[Ilayer];

                        for (let Inode = 0; Inode < numNodesThisLayer; ++Inode) {
                            const thisNode = thisLayer[Inode];

                            thisNode[propertyName] *= scaleFactor;
                        }
                    }
                }
                gatherProperty(propertyName = "") {
                    /**
                     * Returns Array of values for a property across all nodes.
                     * If the property itself happens to be an array, then it
                     pushes each value to the returned array.
                    * Pushes absolute value.
                    * Ignores falsy values.
                    */
                    const net = this.net;
                    const numLayers = this.numLayers;

                    let parDump = new Array();
                    for (let Ilayer = 1; Ilayer < numLayers; ++Ilayer) {
                    
                        const thisLayer = net[Ilayer];
                        const numNodesThisLayer = this.nodesByLayer[Ilayer];
        
                        for (let Inode = 0; Inode < numNodesThisLayer; ++Inode) {
        
                            const thisNode = thisLayer[Inode];
                            parDump.append(thisNode[propertyName]);
                            
                        }
                    }

                    return parDump;
                }
                gatherPropertyAbs_ignore0(propertyName = "") {
                    /**
                     * Returns Array of values for a property across all nodes.
                     * If the property itself happens to be an array, then it
                     pushes each value to the returned array.
                    * Pushes absolute value.
                    * Ignores falsy values.
                    */
                    const net = this.net;
                    const numLayers = this.numLayers;

                    let parDump = new Array();
                    for (let Ilayer = 1; Ilayer < numLayers; ++Ilayer) {
                    
                        const thisLayer = net[Ilayer];
                        const numNodesThisLayer = this.nodesByLayer[Ilayer];
        
                        for (let Inode = 0; Inode < numNodesThisLayer; ++Inode) {
        
                            const thisNode = thisLayer[Inode];
                            parDump.appendAbs_ignore0(thisNode[propertyName]);
                            
                        }
                    }

                    return parDump;
                }

            
            // GENETIC MUTATION FUNCTIONS
                noiseWeights_uniform(span = 0.01) {
                    this.forEachNode( (thisNode) => thisNode.noiseWeights_uniform(span) );
                }
                noiseWeights_sin(span = 0.01) {
                    this.forEachNode( (thisNode) => thisNode.noiseWeights_sin(span) );
                }
                noiseProperty_uniform(span = 0.01, propertyName = "") {
                    this.forEachNode( (thisNode) => thisNode.noiseProperty_uniform(span, propertyName) );
                }
                noiseProperty_sin(span = 0.01, propertyName = "") {
                    this.forEachNode( (thisNode) => thisNode.noiseProperty_sin(span, propertyName) );
                }

                mitose() {
                    /**
                     * Copies the F parameters from `this`
                     * to returned `daughter` NN with no mutation.
                    **/
                    const net = this.net;
                    
                    const daughter = new NN(this.nodesByLayer, 0);

                    daughter.forEachNode( (thisNode, Ilayer, Inode) => {
                        const parentNode = net[Ilayer][Inode];
                        
                        thisNode.weights = parentNode.weights.copy();
                        thisNode.bias    = parentNode.bias;
                        thisNode.α       = parentNode.α;
                    });
                    
                    return daughter;
                }

        }
        
            class Node {
                constructor(Ilayer = 0, Inode = 0, numInputs = 0, weights = new Float32Array(0), bias = 0, toInitialize = true) {
                    // For debugging and ID
                    this.Ilayer     = Ilayer;
                    this.Inode      = Inode;

                    // Number of inputs that this node listens to.
                    this.numInputs  = numInputs;

                    // toInitialize can be set to falsy value if we want to create a new Node
                        // for which we already know the values without spending time on the rand() function
                    // Normalize weights so outputs don't blow up.
                    this.weights    = (weights.length == numInputs) ? weights
                                    : toInitialize ? Float32Array.from (
                                        {length: numInputs},
                                        () => ( (2 * rand() - 1) / this.numInputs)
                                    )
                                    : new Float32Array(numInputs);
                    // Init bias near 0 to give every neuron a fighting chance
                    this.bias       = bias ? bias
                                    : toInitialize ? rand() * 0.2 - 0.1
                                    : 0;
                    // Parameter for leaky activation function
                    this.α          = toInitialize ? rand() * 0.1 - 0.5
                                    : 0;
                    
                    // Forward Prop
                    this.z          = 0;                                    // Weighted sum of inputs, = inputs dotproduct weights
                    this.actIn      = 0;                                    // Input to activation function, = z + bias
                    this.out        = 0;                                    // Output to next Layer
                    
                    // Backprop
                    this.dC_dA      = 0;                                    // D of total cost function w respect to output, = sum of dC_dI s from next layer.
                    this.dA_dB      = 0;                                    // Derivative of output with respect to bias
                    this.dC_dB      = this.dC_dA * this.dA_dB;

                    this.dA_dα      = 0;    	                            // D of A with respect to leaky ReLU parameter
                    this.dC_dα      = 0;
                    
                    this.dC_dZ      = this.dC_dB;                           // Derivative of cost with respect to weighted sum of inputs
                    this.dZ_dW      = new Float32Array(numInputs).fill(0);  // D of weighted sum with respect to each weight
                    this.dC_dW      = new Float32Array(numInputs).fill(0);  // = dC_dZ * dZ_dW
                    this.dC_dI      = new Float32Array(numInputs).fill(0);  // D of cost with respect to prev layer's outputs
                    
                    // For pruning
                    this.Ᾱ  = 0;
                    this.ΔC = 0;
                }
                
                // REINITIALIZATION
                // Ensures that output is between [0, 1] if all inputs are between [0, 1]
                    initWeights_uniform() {
                        this.weights = Array.from(
                            {length: this.numInputs},
                            (x) => ( (2 * rand() - 1) / this.numInputs )
                        );
                    }
                    // Put bias in [-0.1, 0.1] 
                    randBias() {
                        this.bias = rand() * 0.2 - 0.1;
                    }

                // GENETIC TRAINING
                    noiseWeights_uniform(span = 0.01) {
                        this.weights = Float32Array.from(
                            this.weights,
                            (weight) => weight + span * (rand() - 0.5)
                        );
                    }
                    noiseWeights_sin(span = 0.01) {
                        this.weights = Float32Array.from(
                            this.weights,
                            (weight) => weight + span/2 * sin(τ*rand())**3
                        );
                    }
                    noiseProperty_uniform(span = 0.01, propertyName = "") {
                        this[propertyName] += span * (rand() - 0.5);
                    }
                    noiseProperty_sin(span = 0.01, propertyName = "") {
                        this[propertyName] += span/2 * sin(τ*rand())**3;
                    }
                
                // USAGE
                    // Forward Prop
                        // Input layer is passed by reference to save compute
                    F(inputs = new Float32Array()) {
                        this.z      = this.weights.dot(inputs);
                        this.actIn  = this.z + this.bias;
                        this.out    = this.A(this.actIn);
                    }
                    // TRAINING WITH BACKPROP
                    // Back prop
                        // This node is responsible for adding to preceding layer's dC_dA
                    G(inputs = new Float32Array()) {
                    //  this.dC_dA is incremeted by succeeding layer
                        this.dA_dB = this.A_prime(this.actIn);
                        this.dC_dB = this.dC_dA * this.dA_dB;
                        
                        this.dA_dα = this.calc_dA_dα(this.actIn);
                        this.dC_dα = this.dC_dA * this.dA_dα;

                        this.dC_dZ = this.dC_dB;
                        this.dZ_dW = Float32Array.from(inputs);
                        this.dC_dW = this.dZ_dW.scale(this.dC_dZ);

                        // We don't have control over the input layer LOL
                        if (this.Ilayer > 1) {
                        this.dC_dI = this.weights.scale(this.dC_dZ);
                        return this.dC_dI;
                        }
                    }
                    // Gradient Descent
                    D() {
                        this.bias -= this.dC_dB;
                        this.α    -= this.dC_dα;
                        this.weights.minusEquals_ignore0(this.dC_dW);   // Don't step down pruned weights.
                    }

                // Activation Function - Parametric Leaky ReLU with quadratic smoothener from [-0.1 to 0.1]
                    A(x = 0) {
                        if (x >  0.1) return x;
                        if (x < -0.1) return this.α * x;
                        
                        let α = this.α;

                        let a = 2.5   * (1 - α);
                        let b = 0.5   * (1 + α);
                        let c = 0.025 * (1 - α);

                        return a * (x * x) + b * x + c;
                    }
                    A_prime(x = 0) {
                        if (x >  0.1) return 1;
                        if (x < -0.1) return this.α;
                        
                        let α = this.α;

                        let a2 = 5.0 * (1 - α);
                        let b  = 0.5 * (1 + α);

                        return a2 * x + b;
                    }
                    calc_dA_dα(x = 0) {
                        if (x >  0.1) return 0;
                        if (x < -0.1) return x;
                        
                        return -2.5 * (x * x) + 0.5 * x - 0.025;
                    }
            }

/** EXECUTION **/

let myGame = new Game();

/** BROWSWER GAME **/
{
    // myGame.spawn();
    // myGame.spawn();
    // let displayPRE = document.querySelector("#GameState");
    //     displayPRE.innerText = myGame.getGridString();

    // window.onkeydown = function(ev) {
    //     if (!ev.key.includes("Arrow")) return;
    //     log(ev.key);

    //     let gridBefore = Int8Array.from(myGame.grid);
    //     switch (ev.key) {
    //         case "ArrowUp"   : myGame.combineUp();    break;
    //         case "ArrowDown" : myGame.combineDown();  break;
    //         case "ArrowLeft" : myGame.combineLeft();  break;
    //         case "ArrowRight": myGame.combineRight(); break;
    //     }

    //     if (!myGame.grid.equals(gridBefore)) {
    //         myGame.spawn();
    //     }

    //     displayPRE.innerText = myGame.getGridString();
    // }
}

/** MACHINE LEARNING **/

    // 8. Genetic Selection
    {
        const myFam = new Family(undefined, undefined, undefined, undefined);
        const numGens = Infinity;
        
        while (myFam.Igen < numGens) { myFam.stepGen(); }
    }

    // 7. Mitosis & Mutation
    {
        // const motherNN   = new NN([16, 32, 4], 0);
        // const daughters  = new Array();

        // const numDaughters = 5;

        // NNlog("");
        //     NNlog("Mommy Network: ");
        //     NNlog("    W: ", motherNN.gatherProperty("weights").sort( (a, b) => a - b ) );
        //     NNlog("    B: ", motherNN.gatherProperty("bias").sort( (a, b) => a - b ) );
        //     NNlog("    α: ", motherNN.gatherProperty("α").sort( (a, b) => a - b ) );

        // for
        // (
        //     let Idaughter = 0,              span  = 0.1;
        //         Idaughter < numDaughters && span >= 0.01;
        //       ++Idaughter,                  span *= 0.9
        // )
        // {
        //     NNlog("");
        //     NNlog("Daughter", Idaughter, "mitosing with mutation span", span);
        //     const thisDaughter = motherNN.mitose();
        //           thisDaughter.noiseWeights_sin(span);
        //           thisDaughter.noiseProperty_sin(span, "bias");
        //           thisDaughter.noiseProperty_sin(span, "α");
        //           daughters.push(thisDaughter);
        //     NNlog("    W: ", thisDaughter.gatherProperty("weights").sort( (a, b) => a - b ) );
        //     NNlog("    B: ", thisDaughter.gatherProperty("bias").sort( (a, b) => a - b ) );
        //     NNlog("    α: ", thisDaughter.gatherProperty("α").sort( (a, b) => a - b ) );
        // }
    }

    // 6. Letting NN play game
    {
        // let gamerGal = new NNGamer();
        //     gamerLog(gamerGal);
        // while (!gamerGal.makeMove()) {}
    }

    // 5. Picking an output stochastically
    {
        // let distribution = Float32Array.from([0, 0, 0]);
        // let setProbs = Float32Array.from([0.25, 0.25, 0.5]);
        // for (let c = 0; c < 1e5; ++c) {
        //     ++distribution[setProbs.randProbs()];
        // }
        // log(distribution);
    }

    // 4. Updating array elements using map function forEach
    {
        // let myNode = new Node(0, 0, 10);
        // log(myNode.weights);
        // myNode.initWeights_uniform();
        // log(myNode.weights);
    }

    // 3. NN training on my handwritten numerals dataset + Pruning
    
        import * as fs from "fs";
        
        // This code is SOOO TRASH please dont read it
        {
            // function scrambleTraining(trainingData) {
            //     var indexesNotPicked = [];
            //     for (var Ndatum = 0; Ndatum < trainingData.length; Ndatum += 2)
            //         indexesNotPicked.push(Ndatum);
            //     var scrambledTrainingData = [];
            //         for (var Ndatum = 0; Ndatum < trainingData.length; Ndatum += 2) {
            //             var indexOf_indexesNotPicked = Math.floor(Math.random() * indexesNotPicked.length);
            //             var indexToPick = indexesNotPicked[indexOf_indexesNotPicked];
            //             scrambledTrainingData.push(trainingData[indexToPick    ]);
            //             scrambledTrainingData.push(trainingData[indexToPick + 1]);
            //             indexesNotPicked.splice(indexOf_indexesNotPicked, 1);
            //         }
            //     trainingData = scrambledTrainingData;
            //     scrambledTrainingData = undefined; // save RAM
            // }

            // function refresh() {
            //     scrambleTraining(trainData);
            //     inputs          = new Array();
            //     outputs         = new Array();
            //     targetNumerals  = new Array();

            //     for (let Idatum = 0; Idatum < trainData.length; Idatum += 2) {
            //         let digit = trainData[Idatum];
            //         let thisInputVect = Float32Array.from(trainData[Idatum + 1]);
            //         let thisOutputVect = new Float32Array(10);
            //         thisOutputVect[digit] = 1;
                    
            //         inputs.push(thisInputVect);
            //         outputs.push(thisOutputVect);
            //         targetNumerals.push(digit);
            //     }
            // }
            
            // function testAccuracy(netOut, targetNumeral) {
            //     let maxInd = 0;
            //     let maxVal = -Infinity;
            //     netOut.forEach( (v, i) => {
            //         if (v > maxVal) {
            //             maxVal = v;
            //             maxInd = i;
            //         }
            //     });
            //     return (maxInd == targetNumeral) ? 1 : 0;
            // }

            // const tempFileName = "./2048 Master/NN_test.json";
            // let trainData      = JSON.parse(fs.readFileSync(tempFileName).toString().replaceAll("\n", " "));
                
            // let inputs         = new Array();
            // let outputs        = new Array();
            // let targetNumerals = new Array();
            
            // let myNN = new NN([256, 16, 10]);
            //     myNN.learnRate = 0.051;

            // let Iepoch;
            // let sampleCost = Infinity;
            // let accuracy   = 0;
            // let tick = 0, tock = 0;
            // let num_runs = 5;

            // for (Iepoch = -1; Iepoch < num_runs; ++Iepoch) {
            //     refresh();

            //     if (Iepoch >= 0)
            //         NNlog("Epochs Elapsed:", Iepoch, "(", (tock - tick).toFixed(3), "ms )", "Cost", sampleCost.toFixed(5), "Accuracy", (accuracy * 100).toFixed(4), "%");
                
            //     tick = performance.now();
            //     sampleCost = 0;
            //     accuracy = 0;
            //     inputs.forEach( (thisInputVect, Idatum) => {
            //         const thisOutputVect = outputs[Idatum];
            //         myNN.F(thisInputVect);
            //         myNN.G(thisOutputVect);
            //         sampleCost += myNN.cost;
            //         accuracy   += testAccuracy(myNN.netOut, targetNumerals[Idatum]);
            //         if (Iepoch >= 0) {
            //             myNN.D();
            //         }

            //         // NNlog(`At Iepoch ${Iepoch.toString().padStart(3)}, Idatum ${Idatum.toString().padStart(4)}, cost is ${myNN.cost.toPrecision(5).padStart(12)}, output =`,
            //         //     Array.from(myNN.netOut, (x) => x.toPrecision(3).toString().padEnd(12) ).join("")
            //         // );
            //     } );
            //     sampleCost /= inputs.length;
            //     accuracy   /= inputs.length;
            //     // myNN.learnRate += 0.01;
            //     // if (myNN.learnRate > 0.2) myNN.learnRate = 0.2;
            //     tock = performance.now();
                
            // }
            // NNlog("Epochs Elapsed:", Iepoch, "(", (tock - tick).toFixed(3), "ms )", "Cost", sampleCost.toFixed(5), "Accuracy", (accuracy * 100).toFixed(4), "%");
        }
        
        // The same code as above, but with classes so it looks neater :)
        {

            // const tempFileName   = "./2048 Master/NN_test.json";
            // const trainData      = JSON.parse(fs.readFileSync(tempFileName).toString().replaceAll("\n", " "));

            // const inputs         = trainData.sliceStep(2, 1).map( (outVect) => Float32Array.from(outVect) );
            // const outputNumerals = Int8Array.from(trainData.sliceStep(2));
            // const outputs        = new Array();

            // const nodesByLayer   = [256, 32 , 10];
            // const myNN           = new NN(nodesByLayer);
            // const GD             = new Teacher(myNN, inputs, outputs, outputNumerals, 0.05);

            // const myGardener     = new Gardener(myNN, inputs, outputs);

            // GD.setVectorOutputs((outputToken) => {
            //     let outputVect = new Float32Array(10);
            //         outputVect[outputToken] = 1;
            //     return outputVect;
            // });
            
            // // PRUNING
            // const numPruneCycles = Infinity;
            // for (let IpruneCycle = 0; IpruneCycle < numPruneCycles; ++IpruneCycle) {

            //     NNlog("");
            //     NNlog("Prune Cycles Elapsed: ", IpruneCycle);
            //     GD.SGD(Infinity, 0.05, 0.95);

            //     const pruneRatio = 0.1;

            //     let [cost0, acc0] = GD.getCostAccuracy();
            //     NNlog("Before Pruning", pruneRatio*100, "%:");
            //     NNlog("Cost", cost0.toFixed(5), "Accuracy", (acc0 * 100).toFixed(4), "%")
            //     // NNlog("Weights", myNN.gatherProperty("weights"));
                
            //     // myGardener.pruneWeights_magnitude(pruneRatio);   // PASSED
            //     // myGardener.pruneNodes_Ᾱ(pruneRatio);             // PASSED
            //     // myGardener.pruneNodes_ΔC(pruneRatio);            // PASSED
                
            //     let [cost1, acc1] = GD.getCostAccuracy();
            //     NNlog("After Pruning", pruneRatio*100, "%:");
            //     NNlog("Cost", cost1.toFixed(5), "Accuracy", (acc1 * 100).toFixed(4), "%")
            //     // NNlog("Weights", myNN.gatherProperty("weights"));
            // }

        }

    // 2. NN on ONE SAMPLE DATUM
    {
        // let myNN = new NN([16, 8, 4]);
        // let inputLayer = new Float32Array(16).fill(1);
        // let desiredOutput = Float32Array.from([1, 5, 6, 9]);

        
        // // *** NNlog(JSON.stringify(myNN.net));
        // myNN.learnRate = 0.01;
        // let Niter = 0;
        // do {
        //     myNN.F(inputLayer);     // Forward Prop
        //     myNN.G(desiredOutput);  // Back Prop
        //     myNN.D();               // Gradient Descent

        //     // *** NNlog(`At iteration ${Niter.toString().padStart(3)}, cost is ${myNN.cost.toPrecision(5).padStart(12)}, output =`,
        //     // ***     Array.from(myNN.netOut, (x) => x.toPrecision(5).toString().padEnd(12) ).join("")
        //     // *** );
        // } while (++Niter <= 1e5 && myNN.cost > 0.0001);
        // // *** NNlog(JSON.stringify(myNN.net));
    }

    // 1. Node Creation
    {
        // let myNode = new Node(2, 7, 10);
        // log(JSON.stringify(myNode));
    }

/** GAME **/

    // 3 getGridString
    {
        // myGame.grid = Int8Array.from( {length: 16}, () => { return 3 * rand(); });
        // log(myGame.getGridString());
        // myGame.combineUp();
        // log(myGame.getGridString());
    }

    // 2 combining
    {
        // myGame.grid = Int8Array.from( {length: 16}, () => { return 3 * rand(); });
        // myGame.printGrid();
        // myGame.combineUp();
        // myGame.printGrid();
        // myGame.combineDown();
        // myGame.printGrid();
        // myGame.combineLeft();
        // myGame.printGrid();
        // myGame.combineRight();
        // myGame.printGrid();
    }

    // 1) printGrid
    {
        // myGame.grid = Int8Array.from( {length: 16}, () => { return 12 * rand()**4; } );
        // myGame.printGrid();
    }

log("dun!");

// WRITING BREAKS
/***
 *  
 *  Cool Touch
 * 
 *  The night is warm but I'd rather fall into
 *  The cool touch of your arms.
 *  The air desperately tries to take you
 *  dancing through your draping sleeves
 *  but you refuse to leave.
 *  Flowers fall asleep in your hands,
 *  red petals resting over the
 *  bracelets and pens
 *  on your blank skirt.
 *  You're twinkling in the sky
 *  I'm looking for your signs
 *  I'm so scared tonight
 *  I'm meeting you for the first time.
 *  I wish I knew you were shattered to the core
 *  I wish you knew that I loved you for
 *  Now i'm shattered kneeling before
 *  That crescent smile, the lonely mile.
 *  the fact that i got to hold a star for a while.
 *  You're just a star-struck child.
 *  You're just staying for a while.
 *  
 *  Comet
 *  
 *  No, you never needed another.
 *  No, it'll be okay, sister.
 *  Think you can
 *  walk the road alone,
 *  searing cold in your bones
 *  from all the warmth you loaned.
 *  A sky of many, you were the few
 *  who were made for more than
 *  what they have to lose.
 *  But what happens when
 *  the glittery cone is fading
 *  Losing yourself in
 *  what you're fighting far
 *  Give it one more step, one more step
 *  turns to dust just like the rest
 *  please don't take one more step
 *  please don't be like the
 *  
 *  
***/