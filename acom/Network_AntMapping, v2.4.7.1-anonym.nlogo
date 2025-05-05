breed [nodes node]
breed [workers worker]
breed [solutions solution]

workers-own [ memory carry is-carrying? is-active?]
nodes-own [ vector avg-local-distances is-empty?
            ; For GA
            avg-local-distances-GA
          ]
solutions-own [mapping vectors fitness]

globals [ color-list carry-hops number-of-food
          start-modularity start-same-color-links-percentage start-avg-percent-same-color-neighbor t-to-conv start-corrected-avg-percent-same-color-neighbor
          start-avg-distance start-corr-avg-distance
          ; Further recordings
          wall-clock score-before-offload
          ; For GA
          original-mapping original-empty-state best-score iter
          ; For re-run administration
          experiment-num
          ; For checking solution
          problem-state solution-state
          ; For problem analysis
          vector-distances lower-estimate
        ]

to setup
  clear-all

  ask patches [ set pcolor white ]

  ifelse (Topology = "Spatially Clustered") [
    setup-nodes
    setup-spatially-clustered-network
  ]
  [ ifelse (Topology = "Grid-4") [
      set number-of-nodes (int sqrt number-of-nodes) ^ 2
      setup-nodes
      setup-grid4-network
    ]
    [ ifelse (Topology = "Grid-8") [
        setup-nodes
        setup-grid8-network 1
      ]
      [ ifelse (Topology = "Random") [
        setup-nodes
        setup-random-network
        ]
        [ ifelse (Topology = "Watts-Strogatz-1D") [
          setup-nodes
          setup-1DWS-network
          ]
          [
            ifelse (Topology = "Watts-Strogatz-2D") [
              setup-nodes
              setup-2DWS-network
            ]
            [ ifelse (Topology = "Barabási-Albert") [
              setup-nodes
              setup-BA-network
              ]
              [
                print "Illegal network type selected. Using random."
                setup-nodes
                setup-random-network
              ]
            ]
          ]
        ]
      ]
    ]
  ]


  setup-workers

  set number-of-food number-of-nodes * food-density

  set color-list sublist base-colors 0 number-of-colors ; preserved only for compatibility with old modularity calculations
  ask n-of number-of-food nodes [
;    set color one-of color-list
    set vector n-values vdim [random vsize]
    ; TODO: Better color representation
    ; TODO: Make it work with lower dimensions as well
    set color (map-node-color vector)
    set is-empty? false
  ]

  ; For CACHING Purposes
  ask nodes [ calc-avg-local-distances ]

  save-the-problem

  set carry-hops 0
  set start-modularity modularity
  set start-same-color-links-percentage (num-same-color-links) / count links
  set start-avg-percent-same-color-neighbor report-avg-percent-same-color-neighbor
  set start-corrected-avg-percent-same-color-neighbor report-avg-corrected-same-color-neighbor

  set start-avg-distance avg-distance
  set start-corr-avg-distance corr-avg-distance

  set t-to-conv -1
  reset-ticks

  ; If it is a behaviorspace experiment, we save the world (so that we can replicate things later)
  if (length behaviorspace-experiment-name > 0) [ save-the-world ]
end

to-report dist [a b]
  let l min list (length a) (length b)
  let x (map -  a b)
  set x map [i -> i * i] x
  report sqrt sum x
end

to-report create-null-vector
  ifelse  (null-vector-in-middle) [
    report n-values vdim [ vsize / 2 ]
  ][
    report n-values vdim [0]
  ]
end

to-report map-node-color [v]
  report approximate-rgb (item 0 v) (item 1 v) (item 2 v)
end

to setup-nodes
  set-default-shape nodes "circle"
  create-nodes number-of-nodes
  [
    set size 1.5
    set color black
    set vector create-null-vector
    set is-empty? true
  ]
end

to setup-random-network
  let num-links (average-node-degree * number-of-nodes) / 2

  while [count links < num-links] [
    ask one-of nodes [
      let target one-of other nodes

      if (not member? target link-neighbors) [
        create-link-with target
      ]
    ]
  ]

  ; make the network look a little prettier
  repeat 10
  [
    layout-spring nodes links 0.3 (world-width / (sqrt number-of-nodes)) 3
  ]
end

to setup-1DWS-network
  setup-ring-network param-k
  add-shortcuts
end

to setup-2DWS-network
  setup-grid8-network param-k
  add-shortcuts
end

to setup-BA-network
  if (number-of-nodes > average-node-degree) [   ; not an empty network
    ; Creating the core (of average-node-degree)
    let node-list sort nodes
    let core sublist node-list 0 (average-node-degree + 1)
    set node-list sublist node-list (average-node-degree + 1) number-of-nodes

    foreach core [ x ->
      foreach core [ y ->
        if (x != y) [
          ask y [ create-link-with x ]
        ]
      ]
    ]

    ; Creating 'average-node-degree' links for each of the not-core nodes (number-of-nodes - average-node-degree)
    foreach node-list [ x ->
      ask x [
        ; Save a copy of the current set of links. (Simple assignment is not enough as 'links' is special and allowed to grow.
        let old-links link-set links
        while [count link-neighbors < average-node-degree] [
          create-link-with [one-of both-ends] of one-of old-links
        ]
      ]
    ]

    ; make the network look a little prettier
    repeat 10
    [
      layout-radial nodes links first core
    ]
  ]
 end

to setup-spatially-clustered-network
  ask nodes [
    ; for visual reasons, we don't put any nodes *too* close to the edges
    setxy (random-xcor * 0.975) (random-ycor * 0.975)
  ]

  let num-links (average-node-degree * number-of-nodes) / 2
  while [count links < num-links ]
  [
    ask one-of nodes
    [
      let choice (min-one-of (other nodes with [not link-neighbor? myself])
                   [distance myself])
      if choice != nobody [ create-link-with choice ]
    ]
  ]

  ; make the network look a little prettier
  repeat 10
  [
    layout-spring nodes links 0.3 (world-width / (sqrt number-of-nodes)) 1
  ]
end

to setup-grid8-network [ k ]
  let min-id min [who] of nodes
  let side int sqrt number-of-nodes
  let x-padding max-pxcor / side * 2
  let y-padding max-pycor / side * 2

  ask nodes [
    let id who - min-id

    setxy ((id mod side) * x-padding + min-pxcor + 1) ((int (id / side)) * y-padding + min-pycor + 1)
    ;set label (word (id mod side) "," (int (id / side)))
  ]

  ask nodes [
    let x (who - min-id) mod side
    let y int ((who - min-id) / side)

    let mates other nodes with [
      ( (abs( (who - min-id) mod side - x ) <= k) or (abs( (who - min-id) mod side - x ) >= (side - k)) ) and
      ( (abs( int ((who - min-id) / side) - y ) <= k) or (abs( int ((who - min-id) / side) - y) >= (side - k)) )
    ]

    if mates != nobody [ create-links-with mates ]

    ; Wrap around links over the sides of the grid
;    if (x = 0) [
;      set mates other nodes with [ (int ((who - min-id) / side) = y) and ((who - min-id) mod side = (side - 1)) ]
;      if mates != nobody [ create-links-with mates ]
;    ]
;
;    if (y = 0) [
;      set mates other nodes with [ (int ((who - min-id) / side) = (side - 1)) and ((who - min-id) mod side = x) ]
;      if mates != nobody [ create-links-with mates ]
;    ]
  ]
end

to setup-grid4-network
  let min-id min [who] of nodes
  let side int sqrt number-of-nodes
  let x-padding max-pxcor / side * 2
  let y-padding max-pycor / side * 2

  ask nodes [
    let id who - min-id

    setxy ((id mod side) * x-padding + min-pxcor + 1) ((int (id / side)) * y-padding + min-pycor + 1)
    ;set label (word (id mod side) "," (int (id / side)))
  ]

  ask other nodes [
    let x (who - min-id) mod side
    let y int ((who - min-id) / side)

    let mates other nodes with [
      ( (abs( (who - min-id) mod side - x ) < 1) or (abs( (who - min-id) mod side - x ) > (side - 1)) ) xor
      ( (abs( int ((who - min-id) / side) - y ) < 1) or (abs( int ((who - min-id) / side) - y) > (side - 1)) )
    ]

    if mates != nobody [ create-links-with mates ]
  ]
end

to setup-ring-network [ k ]
  ask nodes [
    let id who
    let mates other nodes with [
      (abs( who - id) <= k) or (abs(who - id) >= (number-of-nodes - k))
    ]

    if mates != nobody [ create-links-with mates ]
  ]

  layout-circle sort turtles 15
end

to add-shortcuts
  foreach (range number-of-nodes) [
    x -> foreach (remove x range number-of-nodes) [
      y -> create-stochastic-shortcut x y
    ]
  ]
end

to create-stochastic-shortcut [ x y ]
  if (link x y = nobody) and (random-float 1.0 < param-w / (number-of-nodes * number-of-nodes)) [ ask turtle x [ create-link-with turtle y] ]
end


  to init-worker
    let n one-of nodes
    setxy [xcor] of n  [ycor] of n

    set color white
    set carry create-null-vector
    set memory []
    set is-carrying? false
    set is-active? true
  end


to setup-workers
  set-default-shape workers "bug"
  create-workers number-of-workers
  [
    init-worker
  ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go-as-schelling [ n neighbor-nodes ]

    let local-color [color] of n
    let num-neighbors count neighbor-nodes

    if (color = white) and not (local-color = black) [  ; potential pick-up
      let my-color local-color ; was color
      let num-same-color count neighbor-nodes with [ color = my-color]

      ifelse (behavior = "Threshold") [
        if num-same-color < (threshold * num-neighbors) [
          set color local-color
          ask n [set color black]
        ]
      ] [
        if (random num-neighbors >= num-same-color) [
          set color local-color
          ask n[set color black]
        ]
      ]

    ]

    if not (color = white) and (local-color = black) [ ; potential put-down
      let my-color color
      let num-same-color count neighbor-nodes with [ color = my-color ]

      ifelse (behavior = "Threshold") [
        if (num-same-color >= (threshold * num-neighbors)) [
          ask n [set color my-color]
          set color white
        ]
      ] [
        if (random num-neighbors < num-same-color) [
          ask n[set color my-color]
          set color white
        ]
      ]
    ]
end

to go-as-classic [ n ]
  let local-color [color] of n
  set memory fput local-color memory
  if length memory > memory-length [
    set memory butlast memory
  ]


    if (color = white) and not (local-color = black) [  ; potential pick-up
      let num-same-color length filter [ i -> i = local-color] memory
      let prob k-plus / (k-plus + num-same-color)
      set prob prob * prob

      if (random-float 1 < prob) [
        set color local-color
        ask n [set color black]
      ]
    ]

    if not (color = white) and (local-color = black) [ ; potential put-down
      let my-color color
      let num-same-color length filter [ i -> i = color] memory
      let prob num-same-color / (k-minus + num-same-color)
      set prob prob * prob

      if (random-float 1 < prob) [
        ask n [set color my-color]
        set color white
      ]
    ]
end

to go-as-mapping [ n neighbor-nodes]
  ; We use cache now: ask n [ calc-avg-local-distances ]
  let status-quo [avg-local-distances] of n
  ; to-debug TURN OFF CACHING
  set status-quo calc-distances-from ([vector] of n) neighbor-nodes
  let alternative calc-distances-from carry neighbor-nodes

  let test alternative < status-quo
  if (equal-swaps) [
    set test alternative <= status-quo
  ]

  if (test) [
    let c carry
    let is-c is-carrying?
    let col color
    if (is-c = false) [
      set col black
    ]
    set carry [vector] of n
    set is-carrying? not [is-empty?] of n
    ifelse is-carrying? [
      set color [color] of n
    ][
      set color white
    ]

    ask n [
      ;if is-c [ set vector c ]
      set avg-local-distances alternative ; THIS IS CACHING
      set vector c
      set color col; map-node-color vector
      set is-empty? not is-c
    ]
  ]
end


to go
  ask workers with [is-active? = true] [
    let n one-of nodes-here
    let neighbor-nodes turtle-set n

    let nn neighbor-nodes
    let new-nodes []
    repeat neighbor-distance [
      ask nn [
        let neighs link-neighbors with [not member? self neighbor-nodes]
        set new-nodes (turtle-set new-nodes neighs)
        set neighbor-nodes (turtle-set neighbor-nodes neighs)
      ]
      set nn new-nodes
      set new-nodes []
    ]

    ; removing n from the list
    ask n [
      set neighbor-nodes other neighbor-nodes
    ]

    ; DEBUG
    ;print n
    ;print neighbor-nodes
    ;print [who] of      neighbor-nodes
    ;ask n [set color blue]
    ;ask neighbor-nodes [set color red]

    ifelse (model-version = "Schelling") [
      go-as-schelling n neighbor-nodes
    ] [
      ifelse (model-version = "Classic") [
        go-as-classic n
      ] [
        ifelse (model-version = "Mapping") [
          go-as-mapping n neighbor-nodes
        ][
          print "Unknown model version"
        ]
      ]
    ]



    ; Move to a random neighbor (if any)
    let next-node one-of neighbor-nodes
    if not (next-node = NOBODY) [
      move-to next-node
    ]
  ]

;  let carrying count workers with [not (color = white)]
  let carrying count workers with [is-carrying? = true]
  set carry-hops carry-hops + carrying

  if (carrying = 0) and (t-to-conv = -1) [
    set t-to-conv ticks
  ]
  if (carrying > 1) [
    set t-to-conv -1
  ]

  tick
end

; This method is to force the ants to offload their carry at the end of the run.
; This is in order to avoid distorting the problem (i.e., removing items from the world).
; The offloading process WILL NOT attempt to optimize anymore.
; We allow the ants to JUMP
to offload
  ifelse offload-method = "Rapid" [
    offload-rapidly
  ][
    ifelse offload-method = "Gentle" [
      offload-gently
    ][
      print "UNKNOWN Offload method."
    ]
  ]
  update-plots

  ; We set all agents back to active to allow for possible continuation
  ask workers [ set is-active? true ]
end

to forced-putdown
  let target one-of nodes with [is-empty?]
  if (not (target = Nobody)) [
    move-to target
    let c carry
    let col color
    ask target [
      set vector c
      set is-empty? false
      set color col
      calc-avg-local-distances ; THIS IS NEEDED as plotting is based on cached values
    ]
    set is-carrying? false
    set color white
  ]
end

to offload-rapidly
  ask workers with [is-carrying? = true] [
    forced-putdown
  ]
end

to offload-gently
  let ant one-of workers with [is-carrying? = true]
  while [not (ant = Nobody)] [
    ask ant [
      forced-putdown
      set is-active? false
    ]

    set ant one-of workers with [is-carrying? = true]

    if not (ant = Nobody) [
      repeat gentle-period [ go ]

      ; Bugfix in v2.4.5.2
      if ([is-carrying?] of ant = false) [
        set ant one-of workers with [is-carrying? = true]
      ]
    ]
  ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to run-random-swapper
  set original-mapping map [i -> [vector] of i] (sort-on [who] nodes)
  set original-empty-state  map [i -> [is-empty?] of i] (sort-on [who] nodes)

  ; do some plotting
  set-current-plot "GA: Best Score"
  clear-plot
  set-current-plot-pen "default"

  set iter 0
  repeat num-GA-iterations [
    let a one-of nodes
    let b one-of nodes

    let a-alternative calc-distances-from ([vector] of a)  ([link-neighbors] of b)
    let b-alternative calc-distances-from ([vector] of b)  ([link-neighbors] of a)

    if (([avg-local-distances] of a) - a-alternative) + (([avg-local-distances] of b) - b-alternative) > 0 [ ; TODO:, NOTE: This could be wrong if 'a' and 'b' are neighbors
      let v [vector] of b
      let c [color] of b
      let ie [is-empty?] of b

      ask b [
        set vector [vector] of a
        ; set avg-local-distances b-alternative ; -- CACHING does NOT work, as 'a' and 'b' may be neighbors.
        set color [color] of a
        set is-empty? [is-empty?] of a
      ]
      ask a [
        set vector v
        ; set avg-local-distances a-alternative ; -- CACHING does NOT work, as 'a' and 'b' may be neighbors.
        set color c
        set is-empty? ie
      ]
      ask a [calc-avg-local-distances]
      ask b [calc-avg-local-distances]
    ]
    set iter iter + 1


    ; do some plotting
    if GA-plot? [
      plot avg-distance
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-the-run
  reset-timer
  repeat num-iterations-at-eval [ go ]
  set score-before-offload avg-distance
  offload
  set wall-clock timer
end

to evaluate
  run-SSGA
  do-the-run
end

to run-SSGA-with-1m
  set num-GA-iterations 1000000
  run-SSGA
end

to re-evaluate-with-SSGA-1m
  let file user-new-file
  if (is-string? file) [

    set file word file ".SSGA1m"

    if (file-exists? file) [
      file-delete file
    ]

    file-open file

    set experiment-name "ANTvsSSGA experiment, d=0.95" ; We need this as the 'parameter' will always be reloaded
    foreach (range 1 experiment-num-of-runs) [
      x -> load-the-world experiment-name (x)
      set experiment-num x
      run-SSGA-with-1m
      file-print best-score
      file-flush
      set experiment-name "ANTvsSSGA experiment, d=0.95" ; We need this as the 'parameter' will always be reloaded
    ]

    file-close
  ]
end


to re-run-10x
  let file user-new-file
  if (is-string? file) [

    let original file

    set file word file ".more_runs"

    let i 1
    while [file-exists? file] [
      set file (word original "(" i ").more_runs")
        set i (i + 1)
    ]

;    if (file-exists? file) [
;      ;file-delete file
;    ]

    file-open file

    ;set experiment-name "ANTvsSSGA, vsize=256,  experiment across d" ; We need this as the 'parameter' will always be reloaded
    set experiment-name "ANTvsSSGA, vsize=256,  2023 experiment across d, middle=T, equal=F" ; We need this as the 'parameter' will always be reloaded
    set experiment-num-of-runs 210
;    foreach range experiment-num-of-runs [
    foreach (range experiment-num-of-runs) [
      x ->
        foreach (range 10) [
          y -> load-the-world experiment-name (x + 1)
          set experiment-num x * 10 + y + 1
          random-seed y
          ; FOR DEBUG set num-iterations-at-eval 1000
          do-the-run
          file-print (word avg-distance " " wall-clock " " score-before-offload)
          file-flush
          set experiment-name "ANTvsSSGA, vsize=256,  2023 experiment across d, middle=T, equal=F" ; We need this as the 'parameter' will always be reloaded
          set experiment-num-of-runs 210
        ]
    ]

    file-close
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Based on: https://en.wikipedia.org/wiki/Modularity_(networks)
to-report modularity
  let m count links
  let s num-same-color-links

  set s s / (2 * m)

  let k-i-square map [ x -> report-k-i-squares x m ] lput black color-list

  report (s - sum k-i-square)
end


  ; Helper reporter
  to-report report-k-i-squares [x m]
    let sum-degree sum [link-neighbors] of nodes with [color = x]
    report (sum-degree * sum-degree) / (4 * m * m)
  end

to-report num-same-color-links
  let link-ends [both-ends] of links
  let color-ends map [ x -> [color] of x] link-ends

  let s 0
  foreach color-ends [ x ->
    let f first x
    let l last x
;    if (f != black) and (l != black) and (f = l) [
    if (f = l) [
      set s (s + 1)
    ]
  ]
  report s
end

to-report report-avg-percent-same-color-neighbor
  let s 0

  ask nodes [
    let local-color color
    let num-neighbors count link-neighbors
    if (num-neighbors > 0) [
      let num-same-color count link-neighbors with [ color = local-color]
      set s s + (num-same-color / num-neighbors)
    ]
  ]

  report s / count nodes
end

to-report report-avg-corrected-same-color-neighbor
  let score report-avg-percent-same-color-neighbor
  let c count nodes

  ;let not-carrying count workers with [(color = white)]
  let not-carrying count workers with [is-carrying? = false]
  ; let carrying number-of-workers - not-carrying

  set score score * c
  ;set score score + (carrying * 0)
  set score score + not-carrying ; (not-carrying * 1)

  report score / (c + number-of-workers)
end


to-report %carry-hops
  report carry-hops / number-of-food
end

to-report D-modularity
  report modularity - start-modularity
end

to-report D-same-links
  report 100 * ((num-same-color-links / count links) - start-same-color-links-percentage)
end

to-report D-same-neighbors
  report 100 * (report-avg-percent-same-color-neighbor -  start-avg-percent-same-color-neighbor)
end

to-report D-corrected-same-neighbors
  report 100 * (report-avg-corrected-same-color-neighbor -  start-corrected-avg-percent-same-color-neighbor)
end

to-report %carrying
  report 100 * (count workers with [is-carrying?]) / number-of-workers
end

to-report %same-color-links
  report 100 * num-same-color-links / count links
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report calc-distances-from [v n]
  let s 0
  let num-neighbors count n
  if (num-neighbors > 0) [
    set s sum map [ i -> dist v i] [vector] of n
    set s s / num-neighbors
  ]
  report s
end

to calc-avg-local-distances
 set avg-local-distances (calc-distances-from vector link-neighbors)
end

to-report avg-distance
  ; We DID CACHING, so it should not be needed to: ask nodes [ calc-avg-local-distances ]
  let total-of-avg-distances sum [avg-local-distances] of nodes
  report total-of-avg-distances / count nodes
end

to-report corr-avg-distance
  let score avg-distance
  let c count nodes

  let carrying count workers with [(is-carrying? = true)]

  set score score * c
  let max-distance vdim * vsize
  set score score + carrying * max-distance

  report score / (c + number-of-workers)
end

to-report D-avg-distance
  report avg-distance - start-avg-distance
end

to-report D-corr-avg-distance
  report corr-avg-distance - start-corr-avg-distance
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to export-nodes-color
  let file user-new-file

  if (is-string? file) [

    ;set file word file ".csv"

    if (count nodes > 0) [

      if (file-exists? file) [
        file-delete file
      ]

      file-open file

      let min-id min [who] of nodes
      let side int sqrt number-of-nodes

      foreach range side [
        x -> export-row-color side min-id x
        file-print ""
      ]

      file-close
    ]
  ]
end

to export-row-color [side min-id x]
  foreach range side [
    y -> ask node (y * side + x + min-id ) [file-type color file-type " "]
  ]
end

to export-network
  let file user-new-file
  export-network-with-name file
end

to export-network-with-name [file]
  if (is-string? file) [

    set file word file ".net"

    let n count nodes
    if (n > 0) [

      if (file-exists? file) [
        file-delete file
      ]

      file-open file

      file-print n

      ask links  [
        ask both-ends [ file-type who file-type " " ]
        file-print ""
      ]

      file-close
    ]
  ]
end

to export-vectors
  let file user-new-file
  export-vectors-with-name file
end

to export-vectors-with-name [file]

  if (is-string? file) [

    set file word file ".txt"

    if (count nodes > 0) [

      if (file-exists? file) [
        file-delete file
      ]

      file-open file


      foreach (sort-on [who] nodes) [
        x -> file-type [vector] of x
        file-print ""
      ]

      file-close
    ]
  ]
end

to pick-directory
  set directory-to-save user-directory
end

to save-the-world
  export-world (word directory-to-save "/" behaviorspace-experiment-name " - Run #" behaviorspace-run-number ".state" )
end

to load-the-world [ prefix run-number ]
  import-world (word directory-to-save "/" prefix " - Run #" run-number ".state")
end


; changes all vectors to red (except empty locations)
to all-red
  ask nodes with [is-empty? = false] [
    set vector n-values vdim [ 0 ]
    set vector sentence vsize but-first vector
    set color (map-node-color vector)
  ]
  ; We need to do it separately, to ensure all nodes are set up properly
  ask nodes [ calc-avg-local-distances]
  save-the-problem
end

; changes all vectors to 3 colors (red, green and blue), in equal numbers (except empty locations)
to all-3-colors
  ask nodes with [is-empty? = false] [
    set vector n-values vdim [ 0 ]
    let index random 3
    set vector replace-item index vector vsize
    set color (map-node-color vector)
  ]
  ; We need to do it separately, to ensure all nodes are set up properly
  ask nodes [ calc-avg-local-distances]
  save-the-problem
end

to half-ants
  ;needs forced putdown

  let num count workers with [is-active? = true]
  set num num / 2

  ask n-of num workers with [is-active? = true] [
    forced-putdown
    set is-active? false
  ]
end

to reset-workers
  ; We make sure no vectors will disappear
  offload-rapidly
  ask workers [ die ]
  random-seed RNG-seed
  setup-workers
end

to add-workers
  create-workers num-of-new-workers [
    init-worker
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to save-the-problem
  set problem-state [vector] of nodes with [is-empty? = false]
end

to check-solution
  ifelse not ((length problem-state) = (count nodes with [is-empty? = false])) [
    set solution-state (word "Wrong number of empty locations." (length problem-state) (count nodes with [is-empty? = false]))
  ][
    let vects [vector] of nodes with [is-empty? = false]
    let ps problem-state
    while [not empty? vects] [
      let v first vects
      let p position v ps
      ifelse (p = false) [
        set solution-state "Wrong solution. Vector mismatch."
        set ps []
        set vects [] ; quits the loop
      ][
        set ps remove-item p ps
        set vects butfirst vects
      ]
    ]

    ifelse (not empty? ps) [
      set solution-state "Wrong solution. Vectors not placed."
    ][
      set solution-state "Good solution."
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  Implementation of (SS)GA
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to-report calc-distances-from-GA [v n vs]
  let s 0
  let num-neighbors count n
  if (num-neighbors > 0) [
    let whos [who] of n
    let neighbor-vectors map [i -> item i vs] whos
    set s sum map [ i -> dist v i] neighbor-vectors
    set s s / num-neighbors
  ]
  report s
end

to calc-avg-local-distances-GA [vs]
  let local-vector item who vs
  set avg-local-distances-GA (calc-distances-from-GA local-vector link-neighbors vs)
end

to-report avg-distance-GA [vs]
  ask nodes [ calc-avg-local-distances-GA vs]
  let total-of-avg-distances sum [avg-local-distances-GA] of nodes
  report total-of-avg-distances / count nodes
end

;;;;;;;;;;;;;;

to init-solution
  let IDs range length original-mapping
  set mapping []
  repeat length original-mapping [
    let pick one-of IDs
    set mapping lput pick mapping
    set IDs remove pick IDs
  ]

  set vectors []
  foreach mapping [
    i -> set vectors lput (item i original-mapping ) vectors
  ]

  set fitness avg-distance-GA vectors
end

;;;;;;;;;

;to copy-with-crossover [best other-parent]
;  set mapping []
;
;  foreach range (length original-mapping) [
;    x ->
;    ifelse (random-float 1) < p-crossover [
;      set mapping lput (item x best) mapping
;    ][
;      set mapping lput (item x other-parent) mapping
;    ]
;  ]
;end

to copy-with-crossover [best other-parent]
  let co-point random length original-mapping

  set mapping sublist best 0 co-point ; Note this leaves AT LEAST on for other-parent

  let others []
  foreach other-parent [
    x -> ifelse not member? x mapping [ set mapping lput x mapping ] [ set others lput x others]
  ]
end

to do-point-mutation [loc]
  let targets range length original-mapping
  set targets remove loc targets
  let other-loc one-of targets

  ; We now swap the two mappings
  let temp item other-loc mapping
  set mapping replace-item other-loc mapping item loc mapping
  set mapping replace-item loc mapping temp
end


to mutate
  foreach range (length original-mapping) [
    x -> if (random-float 1) < p-mutation [ do-point-mutation x ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to run-SSGA
  set original-mapping map [i -> [vector] of i] (sort-on [who] nodes)
  set original-empty-state  map [i -> [is-empty?] of i] (sort-on [who] nodes)

  ask solutions [die]

  let p-mut p-mutation
  ifelse (p-mutation < 0) [
    set p-mutation -1 * p-mutation / vdim
  ][
    if (p-mutation >= 1) [
      set p-mutation p-mutation / (number-of-nodes * 1000)
    ]
  ]

  reset-timer
  create-solutions (num-GA-solutions - 1) [
    hide-turtle
    init-solution
  ]

  create-solutions 1 [
    hide-turtle
    set mapping range length original-mapping

    set vectors []
    foreach mapping [
      i -> set vectors lput (item i original-mapping ) vectors
    ]

    set fitness avg-distance-GA vectors
  ]

  ; do some plotting
  set-current-plot "GA: Best Score"
  clear-plot
  set-current-plot-pen "default"

  set iter 0
  repeat num-GA-iterations [
    let best min-one-of solutions [fitness]
    set best-score [fitness] of best
    let worst max-one-of solutions [fitness]

    ask worst [
      ; Putting it in this block, we can use 'other' and thus make sure we do not select the worst as other-parent
      let other-parent one-of other solutions
      copy-with-crossover [mapping] of best [mapping] of other-parent
      mutate

      set vectors []
      foreach mapping [
        i -> set vectors lput (item i original-mapping ) vectors
      ]

      set fitness avg-distance-GA vectors

      set iter iter + 1
    ]


    ; do some plotting
    if GA-plot? [ plot best-score ]
  ]

  set wall-clock timer
  ; restore the original parameter value (it is changed for negative values that has a special meaning: dependence on vdim)
  set p-mutation p-mut
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to copy-best
  let best min-one-of solutions [fitness]
  let best-mapping [mapping] of best

  ask nodes [
    set vector item (item who best-mapping) original-mapping
    set is-empty? item (item who best-mapping) original-empty-state
    set color (map-node-color vector)
  ]
  ; We need to do it separately, to ensure all nodes are restored from backup
  ask nodes [
    calc-avg-local-distances
  ]
end

to copy-orig
  ; This is needed to work for re-initialising ant-algorithms as well (not only SSGA)
  offload-rapidly

  ask nodes [
    set vector item who original-mapping
    set is-empty? item who original-empty-state
    set color (map-node-color vector)
  ]
  ; We need to do it separately, to ensure all nodes are restored from backup
  ask nodes [
    calc-avg-local-distances
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to analyse-problem
  set vector-distances []

  let list-of-vectors [vector] of nodes ; with [is-empty? = false]

  ; This looks complicated, but we want to make sure all pair of nodes is only included once
  foreach list-of-vectors [
    x ->
      let pos position x list-of-vectors
      let sub sublist list-of-vectors (pos + 1) length list-of-vectors
      foreach sub [
        y -> set vector-distances lput (dist x y) vector-distances
      ]
  ]

  let num-links count links
  set vector-distances sort vector-distances
  let lowest sublist vector-distances 0 num-links
  set lower-estimate (sum lowest) / (8 * num-links)

  set-current-plot "Histogram of distances"
  clear-plot
  set-current-plot-pen "default"
  set-histogram-num-bars 100
  set-plot-x-range (min vector-distances) (max vector-distances + 1)
  histogram vector-distances

end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
647
448
-1
-1
13.0
1
10
1
1
1
0
1
1
1
-16
16
-16
16
1
1
1
ticks
30.0

BUTTON
11
14
74
47
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
74
14
137
47
NIL
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
137
14
200
47
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
18
179
190
212
average-node-degree
average-node-degree
1
25
8.0
1
1
NIL
HORIZONTAL

SLIDER
18
212
190
245
number-of-nodes
number-of-nodes
9
1024
256.0
1
1
NIL
HORIZONTAL

SLIDER
20
360
192
393
number-of-workers
number-of-workers
1
250
43.0
1
1
NIL
HORIZONTAL

SLIDER
19
392
191
425
food-density
food-density
0
1.0
0.95
0.05
1
NIL
HORIZONTAL

SLIDER
19
425
191
458
number-of-colors
number-of-colors
1
10
10.0
1
1
NIL
HORIZONTAL

CHOOSER
19
458
190
503
behavior
behavior
"Mapping" "Threshold" "Stochastic"
0

SLIDER
17
524
189
557
threshold
threshold
0
1.0
1.0
0.05
1
NIL
HORIZONTAL

CHOOSER
18
245
190
290
topology
topology
"Grid-4" "Grid-8" "Spatially Clustered" "Random" "Watts-Strogatz-1D" "Watts-Strogatz-2D" "Barabási-Albert"
1

PLOT
663
12
863
162
%Carrying
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot %carrying"

PLOT
663
171
863
321
modularity
NIL
NIL
0.0
10.0
-0.5
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot modularity"

PLOT
663
330
863
480
%same-color-links
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot %same-color-links"

MONITOR
213
454
295
499
#carry-hops
carry-hops
0
1
11

PLOT
866
12
1066
162
avg  of %same-color-neighbor
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot 100 * report-avg-percent-same-color-neighbor"

MONITOR
298
454
384
499
%carry-hops
%carry-hops
3
1
11

MONITOR
387
454
470
499
D-modularity
D-modularity
3
1
11

MONITOR
474
454
556
499
D-same-links
D-same-links
3
1
11

MONITOR
560
454
643
499
D-same-neigh
D-same-neighbors
3
1
11

PLOT
213
503
413
653
degree-distribution
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" "let max-degree max [count link-neighbors] of nodes\nplot-pen-reset  ;; erase what we plotted before\nset-plot-x-range 0 (max-degree + 1)  ;; + 1 to make room for the width of the last bar\nhistogram [count link-neighbors] of nodes\n"

SLIDER
18
289
190
322
param-k
param-k
2
4
4.0
1
1
NIL
HORIZONTAL

SLIDER
18
322
190
355
param-w
param-w
1
100
7.0
1
1
NIL
HORIZONTAL

TEXTBOX
91
339
241
357
x1/N^2\n
11
0.0
1

MONITOR
415
503
528
548
Time-to-Converge
t-to-conv
0
1
11

PLOT
869
174
1069
324
avg-corrected-same-c-neigh
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot 100 * report-avg-corrected-same-color-neighbor"

MONITOR
531
503
647
548
D-corr-same-neigh
D-corrected-same-neighbors
17
1
11

BUTTON
1076
15
1187
48
export-picture
export-nodes-color
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
415
547
553
592
model-version
model-version
"Mapping" "Schelling" "Classic"
0

SLIDER
895
372
1067
405
k-plus
k-plus
0
20
1.0
1
1
NIL
HORIZONTAL

SLIDER
895
406
1067
439
k-minus
k-minus
0
20
15.0
1
1
NIL
HORIZONTAL

SLIDER
895
440
1067
473
memory-length
memory-length
0
100
20.0
1
1
NIL
HORIZONTAL

SLIDER
676
501
848
534
vsize
vsize
1
1000
256.0
1
1
NIL
HORIZONTAL

SLIDER
677
532
849
565
vdim
vdim
3
100
50.0
1
1
NIL
HORIZONTAL

PLOT
1089
75
1289
225
avg-distance
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot avg-distance"
"pen-1" 1.0 0 -13345367 true "" "plot best-score"

PLOT
1088
249
1288
399
corr-avg-distance
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot corr-avg-distance"

MONITOR
1291
73
1388
118
NIL
D-avg-distance
3
1
11

MONITOR
1390
74
1514
119
NIL
D-corr-avg-distance
3
1
11

BUTTON
895
485
980
518
NIL
run-SSGA
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
987
474
1159
507
num-GA-solutions
num-GA-solutions
10
300
100.0
1
1
NIL
HORIZONTAL

SLIDER
988
506
1160
539
num-GA-iterations
num-GA-iterations
1000
1000000
20000.0
100
1
NIL
HORIZONTAL

SLIDER
988
571
1160
604
p-crossover
p-crossover
0
1
0.4
0.05
1
NIL
HORIZONTAL

SLIDER
988
539
1160
572
p-mutation
p-mutation
-10
10000
1000.0
0.005
1
NIL
HORIZONTAL

PLOT
1167
443
1367
593
GA: Best Score
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

BUTTON
896
528
983
561
NIL
copy-best
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
899
572
982
605
NIL
copy-orig\n\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
1265
412
1368
445
GA-plot?
GA-plot?
1
1
-1000

MONITOR
1194
400
1266
445
NIL
best-score
2
1
11

BUTTON
125
53
200
86
NIL
offload
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
11
53
91
86
NIL
evaluate
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1296
127
1381
172
NIL
avg-distance
3
1
11

SLIDER
18
563
191
596
num-iterations-at-eval
num-iterations-at-eval
100
100000
5000.0
100
1
NIL
HORIZONTAL

MONITOR
1137
400
1194
445
#iter
iter
0
1
11

CHOOSER
63
89
201
134
offload-method
offload-method
"Rapid" "Gentle"
1

SLIDER
29
133
201
166
gentle-period
gentle-period
1
100
10.0
1
1
NIL
HORIZONTAL

BUTTON
1187
16
1302
50
NIL
export-vectors
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1302
16
1421
50
NIL
export-network
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1373
424
1480
458
NIL
pick-directory
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1370
365
1586
425
directory-to-save
BehaviorStateGeneratedProblems
1
0
String

BUTTON
1373
580
1553
614
NIL
re-evaluate-with-SSGA-1m
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1372
519
1528
579
experiment-num-of-runs
210.0
1
0
Number

INPUTBOX
1373
459
1800
519
experiment-name
ANTvsSSGA, vsize=256,  2023 experiment across d, middle=T, equal=F
1
0
String

MONITOR
1528
519
1586
564
#run
experiment-num
0
1
11

BUTTON
1372
613
1465
647
NIL
re-run-10x
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1575
45
1642
79
NIL
all-red
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1550
12
1643
46
NIL
all-3-colors
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1370
313
1452
347
NIL
half-ants
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1373
246
1483
280
NIL
reset-workers
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1372
279
1474
313
NIL
add-workers\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
1472
279
1645
312
num-of-new-workers
num-of-new-workers
1
100
15.0
1
1
NIL
HORIZONTAL

INPUTBOX
1482
218
1638
278
RNG-Seed
1.0
1
0
Number

MONITOR
1474
317
1541
362
#workers
count workers with [is-active? = true]
0
1
11

BUTTON
1373
213
1477
247
set-rng-seed
random-seed RNG-Seed
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
22
609
128
643
go-for-period
set original-mapping map [i -> [vector] of i] (sort-on [who] nodes)\nset original-empty-state  map [i -> [is-empty?] of i] (sort-on [who] nodes)\n\nrepeat num-iterations-at-eval [go]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
903
616
1052
650
NIL
run-random-swapper
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
1064
620
1231
640
Uses num-GA-iterations
11
0.0
0

BUTTON
667
597
777
631
NIL
check-solution
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
667
629
889
674
NIL
solution-state
0
1
11

BUTTON
1702
202
1826
236
NIL
analyse-problem
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1702
253
1797
298
NIL
lower-estimate
17
1
11

PLOT
1648
305
1848
455
Histogram of distances
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" "histogram vector-distances"

MONITOR
1543
318
1611
363
#carrying
count workers with [is-carrying? = true]
0
1
11

SWITCH
415
592
580
625
null-vector-in-middle
null-vector-in-middle
0
1
-1000

SWITCH
415
625
537
658
equal-swaps
equal-swaps
0
1
-1000

MONITOR
1240
598
1329
643
GA wall-clock 
wall-clock
3
1
11

SLIDER
28
679
201
712
neighbor-distance
neighbor-distance
1
10
1.0
1
1
NIL
HORIZONTAL

@#$#@#$#@
# Information removed to make the code anonymised
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Grid8, N=225" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>%carrying</metric>
    <metric>t-to-conv</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>report-avg-percent-same-color-neighbor</metric>
    <metric>modularity</metric>
    <metric>%same-color-links</metric>
    <enumeratedValueSet variable="number-of-workers">
      <value value="10"/>
      <value value="30"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Threshold&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="0.1"/>
      <value value="0.15"/>
      <value value="0.2"/>
      <value value="0.25"/>
      <value value="0.3"/>
      <value value="0.35"/>
      <value value="0.4"/>
      <value value="0.45"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Grid8, Good Thresholds, N=225" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>%carrying</metric>
    <metric>t-to-conv</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>report-avg-percent-same-color-neighbor</metric>
    <metric>report-avg-corrected-same-color-neighbor</metric>
    <metric>modularity</metric>
    <metric>%same-color-links</metric>
    <enumeratedValueSet variable="number-of-workers">
      <value value="10"/>
      <value value="30"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Threshold&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0"/>
      <value value="0.125"/>
      <value value="0.25"/>
      <value value="0.375"/>
      <value value="0.5"/>
      <value value="0.625"/>
      <value value="0.75"/>
      <value value="0.875"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="0.1"/>
      <value value="0.15"/>
      <value value="0.2"/>
      <value value="0.25"/>
      <value value="0.3"/>
      <value value="0.35"/>
      <value value="0.4"/>
      <value value="0.45"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Random, degree-sweep, Good Thresholds v2, N=225" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>%carrying</metric>
    <metric>t-to-conv</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>report-avg-percent-same-color-neighbor</metric>
    <metric>report-avg-corrected-same-color-neighbor</metric>
    <metric>modularity</metric>
    <metric>%same-color-links</metric>
    <enumeratedValueSet variable="number-of-workers">
      <value value="10"/>
      <value value="30"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Threshold&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="0"/>
      <value value="0.125"/>
      <value value="0.25"/>
      <value value="0.375"/>
      <value value="0.5"/>
      <value value="0.625"/>
      <value value="0.75"/>
      <value value="0.875"/>
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="WS-1D, Good Thresholds v2, N=225" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>%carrying</metric>
    <metric>t-to-conv</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>report-avg-percent-same-color-neighbor</metric>
    <metric>report-avg-corrected-same-color-neighbor</metric>
    <metric>modularity</metric>
    <metric>%same-color-links</metric>
    <enumeratedValueSet variable="number-of-workers">
      <value value="10"/>
      <value value="30"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Threshold&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Watts-Strogatz-1D&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="0"/>
      <value value="0.125"/>
      <value value="0.25"/>
      <value value="0.375"/>
      <value value="0.5"/>
      <value value="0.625"/>
      <value value="0.75"/>
      <value value="0.875"/>
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="BA, Good Thresholds v2, N=225" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>%carrying</metric>
    <metric>t-to-conv</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>report-avg-percent-same-color-neighbor</metric>
    <metric>report-avg-corrected-same-color-neighbor</metric>
    <metric>modularity</metric>
    <metric>%same-color-links</metric>
    <enumeratedValueSet variable="number-of-workers">
      <value value="10"/>
      <value value="30"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Threshold&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Barabási-Albert&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="0"/>
      <value value="0.125"/>
      <value value="0.25"/>
      <value value="0.375"/>
      <value value="0.5"/>
      <value value="0.625"/>
      <value value="0.75"/>
      <value value="0.875"/>
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Grid8, Stochastic, N=225" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>%carrying</metric>
    <metric>t-to-conv</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>report-avg-percent-same-color-neighbor</metric>
    <metric>report-avg-corrected-same-color-neighbor</metric>
    <metric>modularity</metric>
    <metric>%same-color-links</metric>
    <enumeratedValueSet variable="number-of-workers">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Threshold&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Classic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="5"/>
      <value value="10"/>
      <value value="15"/>
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA experiment" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="250"/>
      <value value="400"/>
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="255"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
      <value value="6"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA experiment d=0.5" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="250"/>
      <value value="400"/>
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="255"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
      <value value="6"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA experiment, d=1.0" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <metric>wall-clock</metric>
    <metric>score-before-offload</metric>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="250"/>
      <value value="400"/>
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="255"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
      <value value="6"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA experiment d=0.95" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="250"/>
      <value value="400"/>
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="255"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
      <value value="6"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA, vsize=256,  experiment across d" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="0.95"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
      <value value="6"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA experiment across M" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="255"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA experiment across M, D=6" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="255"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA experiment across M, D=50" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="255"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA experiment d=0.95, p-mutation=1/N, pt2" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="400"/>
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="255"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="-0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
      <value value="6"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA experiment d=0.95, p-mutation=1divN, pt1" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="255"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="-0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
      <value value="6"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA, v=256, experiment across M, D=50" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA, v=256, experiment across M, D=6" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA. v=256, experiment across M, D=3" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=50, d=0.95" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=6, d=0.95" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=3, d=0.95" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=50, d=0.95, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=6, d=0.95, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=3, d=0.95, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA, vsize=256,  2023 experiment across d, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="0.95"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
      <value value="6"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA, v=256, experiment across M, D=50, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA, v=256, experiment across M, D=6, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA. v=256, experiment across M, D=3, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=50, d=0.95, middle=T, equal=F" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=6, d=0.95, middle=T, equal=F" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=3, d=0.95, middle=T, equal=F" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA, vsize=256,  2023 experiment across d, middle=T, equal=F" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="0.95"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
      <value value="6"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=50, d=0.95, p-mut=1#D, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=3, d=0.95, p-mut=1#D, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=6, d=0.95, p-mut=1#D, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=3, d=0.95, p-mut=1#N, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=6, d=0.95, p-mut=1#N, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ANTvsSSGA,  v=256, experiment vs N, D=50, d=0.95, p-mut=1#N, middle=T, equal=T" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>evaluate</go>
    <timeLimit steps="1"/>
    <metric>avg-distance</metric>
    <metric>best-score</metric>
    <metric>carry-hops</metric>
    <metric>%carry-hops</metric>
    <metric>D-modularity</metric>
    <metric>D-same-links</metric>
    <metric>D-same-neighbors</metric>
    <metric>t-to-conv</metric>
    <metric>D-corrected-same-neighbors</metric>
    <metric>D-avg-distance</metric>
    <metric>D-corr-avg-distance</metric>
    <enumeratedValueSet variable="directory-to-save">
      <value value="&quot;BehaviorStateGeneratedProblems&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behavior">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-workers">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-nodes">
      <value value="225"/>
      <value value="256"/>
      <value value="400"/>
      <value value="1024"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="GA-plot?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-colors">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-w">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gentle-period">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vsize">
      <value value="256"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-version">
      <value value="&quot;Mapping&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="average-node-degree">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-iterations-at-eval">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="offload-method">
      <value value="&quot;Gentle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-plus">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-mutation">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="memory-length">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-solutions">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-minus">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="param-k">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-density">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="topology">
      <value value="&quot;Grid-8&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-crossover">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vdim">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-GA-iterations">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-vector-in-middle">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="equal-swaps">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
