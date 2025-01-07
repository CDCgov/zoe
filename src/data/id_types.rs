pub(crate) trait FastaIDs {
    fn get_id_taxon(&self) -> Option<(&str, &str)>;
}

impl<S: AsRef<str> + ?Sized> FastaIDs for S {
    fn get_id_taxon(&self) -> Option<(&str, &str)> {
        let s = self.as_ref();
        if let Some(end_id) = s.find('{')
            && let Some(start_annot) = s.rfind('{')
            && let Some(annot_offset) = s[start_annot..].find('}')
            && annot_offset > 1
        {
            Some((&s[..end_id], &s[start_annot + 1..start_annot + annot_offset]))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn get_id_taxon_test() {
        let case_result = [
            ("id{clade}", Some(("id", "clade"))),
            ("id2{clade1}{clade2}", Some(("id2", "clade2"))),
            ("id3{{clade}}", Some(("id3", "clade"))),
            ("id4{}", None),
            ("id5_clade", None),
            ("", None),
        ];

        for (input, expected) in case_result {
            assert_eq!(expected, input.get_id_taxon());
        }
    }
}
