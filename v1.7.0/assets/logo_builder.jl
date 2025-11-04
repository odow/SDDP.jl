using Luxor

function logo_with_text()
    P = 20
    R = 100
    H = 10
    Drawing(2R + 2P, 2R + H + 2P, joinpath(@__DIR__, "logo.svg"))
    function draw_box(color, angle, w)
        sethue(color)
        rotate(angle)
        box(Point(w * R, 0), R, 3R, :fill)
        rotate(-angle)
        return
    end
    origin()
    setopacity(3 / 4)
    circle(Point(0, -H / 2 - 50), R, :clip)
    draw_box(Luxor.julia_green, 2π / 3, 1 / 3)
    draw_box(Luxor.julia_purple, π / 3, 2 / 5)
    draw_box(Luxor.julia_red, π / 6, 3 / 4)
    clipreset()
    setopacity(1)
    setcolor("black")
    setfont("Arial", 60)
    settext(
        "<b>SDDP.jl</b>",
        Point(0, R + H / 2 + P);
        halign = "center",
        markup = true,
    )
    finish()
    return
end

function logo_without_text()
    P = 10
    R = 100
    H = 10
    Drawing(2R + 2P, 2R + 2P - 50, joinpath(@__DIR__, "logo_without_text.svg"))
    p = origin(Point(110, R + P - 25 + 55 / 2))
    function draw_box(color, angle, w)
        sethue(color)
        rotate(angle)
        box(Point(w * R, 0), R, 3R, :fill)
        rotate(-angle)
        return
    end
    setopacity(3 / 4)
    circle(Point(0, -H / 2 - 50), R, :clip)
    draw_box(Luxor.julia_green, 2π / 3, 1 / 3)
    draw_box(Luxor.julia_purple, π / 3, 2 / 5)
    draw_box(Luxor.julia_red, π / 6, 3 / 4)
    finish()
    return
end

logo_with_text()
logo_without_text()
