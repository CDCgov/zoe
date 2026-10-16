
use crate::data::sam::{SamOptField, SamOptRaw, SamOptRawView, SamOptValue};

#[test]
fn sam_opt_raw_len() {
    assert_eq!(SamOptRaw::new().len(), 0);
    assert_eq!(SamOptRaw::default().len(), 0);
    assert_eq!(SamOptRaw(String::from(" ")).len(), 1);
    assert_eq!(SamOptRaw(String::from("\t")).len(), 2);
    assert_eq!(SamOptRaw(String::from("\t\t")).len(), 3);

    assert_eq!(SamOptRaw(String::from("field")).len(), 1);
    assert_eq!(SamOptRaw(String::from("field1\tfield2")).len(), 2);
    assert_eq!(SamOptRaw(String::from("field1\tfield2\t")).len(), 3);
    assert_eq!(SamOptRaw(String::from("\tfield1\tfield2")).len(), 3);
    assert_eq!(SamOptRaw(String::from("field1\t\tfield2")).len(), 3);

    assert_eq!(SamOptRawView::new().len(), 0);
    assert_eq!(SamOptRawView::default().len(), 0);
    assert_eq!(SamOptRawView(" ").len(), 1);
    assert_eq!(SamOptRawView("\t").len(), 2);
    assert_eq!(SamOptRawView("\t\t").len(), 3);

    assert_eq!(SamOptRawView("field").len(), 1);
    assert_eq!(SamOptRawView("field1\tfield2").len(), 2);
    assert_eq!(SamOptRawView("field1\tfield2\t").len(), 3);
    assert_eq!(SamOptRawView("\tfield1\tfield2").len(), 3);
    assert_eq!(SamOptRawView("field1\t\tfield2").len(), 3);
}

#[test]
fn sam_opt_raw_is_empty() {
    assert!(SamOptRaw::new().is_empty());
    assert!(SamOptRaw::default().is_empty());
    assert!(!SamOptRaw(String::from(" ")).is_empty());
    assert!(!SamOptRaw(String::from("\t")).is_empty());
    assert!(!SamOptRaw(String::from("\t\t")).is_empty());

    assert!(!SamOptRaw(String::from("field")).is_empty());
    assert!(!SamOptRaw(String::from("field1\tfield2")).is_empty());
    assert!(!SamOptRaw(String::from("field1\tfield2\t")).is_empty());
    assert!(!SamOptRaw(String::from("\tfield1\tfield2")).is_empty());
    assert!(!SamOptRaw(String::from("field1\t\tfield2")).is_empty());

    assert!(SamOptRawView::new().is_empty());
    assert!(SamOptRawView::default().is_empty());
    assert!(!SamOptRawView(" ").is_empty());
    assert!(!SamOptRawView("\t").is_empty());
    assert!(!SamOptRawView("\t\t").is_empty());

    assert!(!SamOptRawView("field").is_empty());
    assert!(!SamOptRawView("field1\tfield2").is_empty());
    assert!(!SamOptRawView("field1\tfield2\t").is_empty());
    assert!(!SamOptRawView("\tfield1\tfield2").is_empty());
    assert!(!SamOptRawView("field1\t\tfield2").is_empty());
}

#[test]
#[allow(clippy::too_many_lines)]
fn sam_opt_raw_iter_raw() {
    assert_eq!(SamOptRaw::new().iter_raw().count(), 0);
    assert_eq!(SamOptRaw::default().iter_raw().count(), 0);
    {
        let opts = SamOptRaw(String::from(" "));
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some(" "));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRaw(String::from("\t"));
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRaw(String::from("\t\t"));
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), None);
    }

    {
        let opts = SamOptRaw(String::from("field"));
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some("field"));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRaw(String::from("field1\tfield2"));
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some("field1"));
        assert_eq!(iter.next(), Some("field2"));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRaw(String::from("field1\tfield2\t"));
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some("field1"));
        assert_eq!(iter.next(), Some("field2"));
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRaw(String::from("\tfield1\tfield2"));
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), Some("field1"));
        assert_eq!(iter.next(), Some("field2"));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRaw(String::from("field1\t\tfield2"));
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some("field1"));
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), Some("field2"));
        assert_eq!(iter.next(), None);
    }

    assert_eq!(SamOptRawView::new().iter_raw().count(), 0);
    assert_eq!(SamOptRawView::default().iter_raw().count(), 0);
    {
        let opts = SamOptRawView(" ");
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some(" "));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRawView("\t");
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRawView("\t\t");
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), None);
    }

    {
        let opts = SamOptRawView("field");
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some("field"));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRawView("field1\tfield2");
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some("field1"));
        assert_eq!(iter.next(), Some("field2"));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRawView("field1\tfield2\t");
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some("field1"));
        assert_eq!(iter.next(), Some("field2"));
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRawView("\tfield1\tfield2");
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), Some("field1"));
        assert_eq!(iter.next(), Some("field2"));
        assert_eq!(iter.next(), None);
    }
    {
        let opts = SamOptRawView("field1\t\tfield2");
        let mut iter = opts.iter_raw();
        assert_eq!(iter.next(), Some("field1"));
        assert_eq!(iter.next(), Some(""));
        assert_eq!(iter.next(), Some("field2"));
        assert_eq!(iter.next(), None);
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn sam_opt_raw_push() {
    {
        let mut opt = SamOptRaw::new();
        opt.push("AS", &SamOptValue::Int(10));
        assert_eq!(opt.to_string(), "AS:i:10");
    }
    {
        let mut opt = SamOptRaw::default();
        opt.push("AS", &SamOptValue::Int(10));
        assert_eq!(opt.to_string(), "AS:i:10");
    }
    {
        let mut opt = SamOptRaw(String::from(" "));
        opt.push("AS", &SamOptValue::Int(10));
        assert_eq!(opt.to_string(), " \tAS:i:10");
    }
    {
        let mut opt = SamOptRaw(String::from("\t"));
        opt.push("AS", &SamOptValue::Int(10));
        assert_eq!(opt.to_string(), "\t\tAS:i:10");
    }
    {
        let mut opt = SamOptRaw(String::from("\t\t"));
        opt.push("AS", &SamOptValue::Int(10));
        assert_eq!(opt.to_string(), "\t\t\tAS:i:10");
    }

    {
        let mut opt = SamOptRaw(String::from("field"));
        opt.push("AS", &SamOptValue::Int(10));
        assert_eq!(opt.to_string(), "field\tAS:i:10");
    }
    {
        let mut opt = SamOptRaw(String::from("field1\tfield2"));
        opt.push("AS", &SamOptValue::Int(10));
        assert_eq!(opt.to_string(), "field1\tfield2\tAS:i:10");
    }
    {
        let mut opt = SamOptRaw(String::from("field1\tfield2\t"));
        opt.push("AS", &SamOptValue::Int(10));
        assert_eq!(opt.to_string(), "field1\tfield2\t\tAS:i:10");
    }
    {
        let mut opt = SamOptRaw(String::from("\tfield1\tfield2"));
        opt.push("AS", &SamOptValue::Int(10));
        assert_eq!(opt.to_string(), "\tfield1\tfield2\tAS:i:10");
    }
    {
        let mut opt = SamOptRaw(String::from("field1\t\tfield2"));
        opt.push("AS", &SamOptValue::Int(10));
        assert_eq!(opt.to_string(), "field1\t\tfield2\tAS:i:10");
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn sam_opt_raw_from_iter() {
    assert_eq!(std::iter::empty::<String>().collect::<SamOptRaw>().0, SamOptRaw::new().0);
    assert_eq!(
        std::iter::once(String::from("field")).collect::<SamOptRaw>().to_string(),
        SamOptRaw(String::from("field")).to_string()
    );
    assert_eq!(
        ["field1", "field2"]
            .into_iter()
            .map(String::from)
            .collect::<SamOptRaw>()
            .to_string(),
        SamOptRaw(String::from("field1\tfield2")).to_string()
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn sam_opt_raw_into_iter() {
    assert_eq!(SamOptRaw::new().into_iter().count(), 0);
    assert_eq!(SamOptRaw::default().into_iter().count(), 0);
    {
        let opts = SamOptRaw(String::from(" "));
        let mut iter = opts.into_iter();
        assert!(iter.next().is_some_and(|res| res.is_err()));
        assert!(iter.next().is_none());
    }
    {
        let opts = SamOptRaw(String::from("\t"));
        let mut iter = opts.into_iter();
        assert!(iter.next().is_some_and(|res| res.is_err()));
        assert!(iter.next().is_some_and(|res| res.is_err()));
        assert!(iter.next().is_none());
    }
    {
        let opts = SamOptRaw(String::from("\t\t"));
        let mut iter = opts.into_iter();
        assert!(iter.next().is_some_and(|res| res.is_err()));
        assert!(iter.next().is_some_and(|res| res.is_err()));
        assert!(iter.next().is_some_and(|res| res.is_err()));
        assert!(iter.next().is_none());
    }

    let expected1 = SamOptField {
        tag:   *b"AS",
        value: SamOptValue::Int(10),
    };
    let expected2 = SamOptField {
        tag:   *b"AS",
        value: SamOptValue::Int(20),
    };

    {
        let opts = SamOptRaw(String::from("AS:i:10"));
        let mut iter = opts.into_iter();
        assert_eq!(iter.next().unwrap().unwrap(), expected1);
        assert!(iter.next().is_none());
    }
    {
        let opts = SamOptRaw(String::from("AS:i:10\tAS:i:20"));
        let mut iter = opts.into_iter();
        assert_eq!(iter.next().unwrap().unwrap(), expected1);
        assert_eq!(iter.next().unwrap().unwrap(), expected2);
        assert!(iter.next().is_none());
    }
    {
        let opts = SamOptRaw(String::from("AS:i:10\tAS:i:20\t"));
        let mut iter = opts.into_iter();
        assert_eq!(iter.next().unwrap().unwrap(), expected1);
        assert_eq!(iter.next().unwrap().unwrap(), expected2);
        assert!(iter.next().is_some_and(|res| res.is_err()));
        assert!(iter.next().is_none());
    }
    {
        let opts = SamOptRaw(String::from("\tAS:i:10\tAS:i:20"));
        let mut iter = opts.into_iter();
        assert!(iter.next().is_some_and(|res| res.is_err()));
        assert_eq!(iter.next().unwrap().unwrap(), expected1);
        assert_eq!(iter.next().unwrap().unwrap(), expected2);
        assert!(iter.next().is_none());
    }
    {
        let opts = SamOptRaw(String::from("AS:i:10\t\tAS:i:20"));
        let mut iter = opts.into_iter();
        assert_eq!(iter.next().unwrap().unwrap(), expected1);
        assert!(iter.next().is_some_and(|res| res.is_err()));
        assert_eq!(iter.next().unwrap().unwrap(), expected2);
        assert!(iter.next().is_none());
    }
}
