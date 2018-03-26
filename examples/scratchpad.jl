mutable struct Happy
    intensity::Int
    reason::String
end

mutable struct Sad
    intensity::Int
    reason::String
    howlong::Int
end

Emotion = Union{Happy, Sad}

function emotion_reader(self::Emotion)

    println(typeof(self))
    println("Intensity: ", getfield(self, :intensity) )
    println("Reason: ", getfield(self, :reason) )

end

a_happy_emotion = Happy(10, "great day")

a_sad_emotion = Sad(7, "bad tea", 2)

emotion_reader(a_happy_emotion)

println("\n")

emotion_reader(a_sad_emotion)