test_that("getSequenceAnnotations works", {
    probe_meta = data.frame(
        PROBE_ID=c("A;1","A;2"),
        PROBE_SEQUENCE = c("MSGSASFEGGVFSPYL","SGSASFEGGVFSPYLT"));
    expect_equal(
        getSequenceAnnotations("A_1_2", probe_meta),
        data.frame(
            Overlap.Seq.Length = 15,
            Full.Seq.Start = 1,
            Full.Seq.Stop = 17,
            Full.Seq.Length = 17,
            First.Seq = "MSGSASFEGGVFSPYL",
            Last.Seq = "SGSASFEGGVFSPYLT",
            Overlap.Seq = "SGSASFEGGVFSPYL",
            Full.Seq = "MSGSASFEGGVFSPYLT"
        )
    )
    expect_equal(
        getSequenceAnnotations("A_1_1", probe_meta),
        data.frame(
            Overlap.Seq.Length = 16,
            Full.Seq.Start = 1,
            Full.Seq.Stop = 16,
            Full.Seq.Length = 16,
            First.Seq = "MSGSASFEGGVFSPYL",
            Last.Seq = "MSGSASFEGGVFSPYL",
            Overlap.Seq = "MSGSASFEGGVFSPYL",
            Full.Seq = "MSGSASFEGGVFSPYL"
        )
    )
})
