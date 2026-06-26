use crate::{
    alignment::phmm::{
        CorePhmm, DomainPhmm, GlobalPhmm, LayerParams, LocalPhmm, SemiLocalPhmm,
        indexing::{
            GetCore, GetCoreMut, GetLayer, GetLayerMut, GetMapping, GetModule, GetModuleMut, GetPartsMut, PhmmIndex,
            PhmmIndexRange, PhmmIndexable,
        },
        modules::{DomainModule, LocalModule, SemiLocalModule},
    },
    data::{
        ByteIndexMap,
        views::{AsView, AsViewMut, AssocOwnedType, AssocViewMutType, AssocViewType, ToOwnedData, ToView},
    },
    prelude::{DataView, DataViewMut},
};

/// The corresponding immutable view type for [`GlobalPhmm`].
///
/// See [Views](crate::data#views) for more details.
#[derive(Eq, PartialEq, Debug)]
pub struct GlobalPhmmView<'a, T, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    core:    &'a CorePhmm<T, S>,
}

/// The corresponding immutable view type for [`LocalPhmm`].
///
/// See [Views](crate::data#views) for more details.
#[derive(Eq, PartialEq, Debug)]
pub struct LocalPhmmView<'a, T, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    core:    &'a CorePhmm<T, S>,
    begin:   &'a LocalModule<T, S>,
    end:     &'a LocalModule<T, S>,
}

/// The corresponding immutable view type for [`DomainPhmm`].
///
/// See [Views](crate::data#views) for more details.
#[derive(Eq, PartialEq, Debug)]
pub struct DomainPhmmView<'a, T, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    core:    &'a CorePhmm<T, S>,
    begin:   &'a DomainModule<T, S>,
    end:     &'a DomainModule<T, S>,
}

/// The corresponding immutable view type for [`SemiLocalPhmm`].
///
/// See [Views](crate::data#views) for more details.
#[derive(Eq, PartialEq, Debug)]
pub struct SemiLocalPhmmView<'a, T, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    core:    &'a CorePhmm<T, S>,
    begin:   &'a SemiLocalModule<T>,
    end:     &'a SemiLocalModule<T>,
}

/// The corresponding mutable view type for [`GlobalPhmm`].
///
/// See [Views](crate::data#views) for more details.
pub struct GlobalPhmmViewMut<'a, T, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    core:    &'a mut CorePhmm<T, S>,
}

/// The corresponding mutable view type for [`LocalPhmm`].
///
/// See [Views](crate::data#views) for more details.
pub struct LocalPhmmViewMut<'a, T, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    core:    &'a mut CorePhmm<T, S>,
    begin:   &'a mut LocalModule<T, S>,
    end:     &'a mut LocalModule<T, S>,
}

/// The corresponding mutable view type for [`DomainPhmm`].
///
/// See [Views](crate::data#views) for more details.
pub struct DomainPhmmViewMut<'a, T, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    core:    &'a mut CorePhmm<T, S>,
    begin:   &'a mut DomainModule<T, S>,
    end:     &'a mut DomainModule<T, S>,
}

/// The corresponding mutable view type for [`SemiLocalPhmm`].
///
/// See [Views](crate::data#views) for more details.
pub struct SemiLocalPhmmViewMut<'a, T, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    core:    &'a mut CorePhmm<T, S>,
    begin:   &'a mut SemiLocalModule<T>,
    end:     &'a mut SemiLocalModule<T>,
}

// Explicit impls for Copy and Clone are needed, since the derived versions only
// implement it when `T` implements Copy and Clone

impl<T, const S: usize> Clone for GlobalPhmmView<'_, T, S> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<T, const S: usize> Copy for GlobalPhmmView<'_, T, S> {}

impl<T, const S: usize> Clone for LocalPhmmView<'_, T, S> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<T, const S: usize> Copy for LocalPhmmView<'_, T, S> {}

impl<T, const S: usize> Clone for DomainPhmmView<'_, T, S> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<T, const S: usize> Copy for DomainPhmmView<'_, T, S> {}

impl<T, const S: usize> Clone for SemiLocalPhmmView<'_, T, S> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<T, const S: usize> Copy for SemiLocalPhmmView<'_, T, S> {}

impl<T, const S: usize> AssocOwnedType for GlobalPhmm<T, S> {
    type Owned = GlobalPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for GlobalPhmm<T, S>
where
    T: Clone + 'static,
{
    type View<'a> = GlobalPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for GlobalPhmm<T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = GlobalPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AssocOwnedType for GlobalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type Owned = GlobalPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for GlobalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type View<'a> = GlobalPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for GlobalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = GlobalPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AssocOwnedType for GlobalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type Owned = GlobalPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for GlobalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type View<'a> = GlobalPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for GlobalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = GlobalPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AsView for GlobalPhmm<T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        GlobalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
        }
    }
}

impl<T, const S: usize> AsViewMut for GlobalPhmm<T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view_mut(&mut self) -> Self::ViewMut<'_> {
        GlobalPhmmViewMut {
            mapping: self.mapping(),
            core:    self.core_mut(),
        }
    }
}

impl<T, const S: usize> ToOwnedData for GlobalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        GlobalPhmm {
            mapping: self.mapping,
            core:    self.core.clone(),
        }
    }
}

impl<'a, T, const S: usize> DataView<'a> for GlobalPhmmView<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        GlobalPhmmView {
            mapping: self.mapping,
            core:    self.core,
        }
    }
}

impl<T, const S: usize> AsView for GlobalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        GlobalPhmmView {
            mapping: self.mapping,
            core:    self.core,
        }
    }
}

impl<'a, T, const S: usize> ToView<'a> for GlobalPhmmViewMut<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_view(self) -> Self::View<'a> {
        GlobalPhmmView {
            mapping: self.mapping,
            core:    self.core,
        }
    }
}

impl<T, const S: usize> ToOwnedData for GlobalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        GlobalPhmm {
            mapping: self.mapping,
            core:    self.core.clone(),
        }
    }
}

impl<'a, T, const S: usize> DataViewMut<'a> for GlobalPhmmViewMut<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b, {
        GlobalPhmmViewMut {
            mapping: self.mapping,
            core:    self.core,
        }
    }
}

impl<T, const S: usize> AssocOwnedType for LocalPhmm<T, S> {
    type Owned = LocalPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for LocalPhmm<T, S>
where
    T: Clone + 'static,
{
    type View<'a> = LocalPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for LocalPhmm<T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = LocalPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AssocOwnedType for LocalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type Owned = LocalPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for LocalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type View<'a> = LocalPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for LocalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = LocalPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AssocOwnedType for LocalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type Owned = LocalPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for LocalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type View<'a> = LocalPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for LocalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = LocalPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AsView for LocalPhmm<T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        LocalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   self.begin(),
            end:     self.end(),
        }
    }
}

impl<T, const S: usize> AsViewMut for LocalPhmm<T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view_mut(&mut self) -> Self::ViewMut<'_> {
        let mapping = self.mapping();
        let (core, begin, end) = self.parts_mut();

        LocalPhmmViewMut {
            mapping,
            core,
            begin,
            end,
        }
    }
}

impl<T, const S: usize> ToOwnedData for LocalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        LocalPhmm::new(self.mapping, self.core.clone(), self.begin().clone(), self.end().clone())
    }
}

impl<'a, T, const S: usize> DataView<'a> for LocalPhmmView<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        LocalPhmmView {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin(),
            end:     self.end(),
        }
    }
}

impl<T, const S: usize> AsView for LocalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        LocalPhmmView {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

impl<'a, T, const S: usize> ToView<'a> for LocalPhmmViewMut<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_view(self) -> Self::View<'a> {
        LocalPhmmView {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

impl<T, const S: usize> ToOwnedData for LocalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        LocalPhmm {
            mapping: self.mapping,
            core:    self.core.clone(),
            begin:   self.begin.clone(),
            end:     self.end.clone(),
        }
    }
}

impl<'a, T, const S: usize> DataViewMut<'a> for LocalPhmmViewMut<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b, {
        LocalPhmmViewMut {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

impl<T, const S: usize> AssocOwnedType for SemiLocalPhmm<T, S> {
    type Owned = SemiLocalPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for SemiLocalPhmm<T, S>
where
    T: Clone + 'static,
{
    type View<'a> = SemiLocalPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for SemiLocalPhmm<T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = SemiLocalPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AssocOwnedType for SemiLocalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type Owned = SemiLocalPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for SemiLocalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type View<'a> = SemiLocalPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for SemiLocalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = SemiLocalPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AssocOwnedType for SemiLocalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type Owned = SemiLocalPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for SemiLocalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type View<'a> = SemiLocalPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for SemiLocalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = SemiLocalPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AsView for SemiLocalPhmm<T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        SemiLocalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   self.begin(),
            end:     self.end(),
        }
    }
}

impl<T, const S: usize> AsViewMut for SemiLocalPhmm<T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view_mut(&mut self) -> Self::ViewMut<'_> {
        let mapping = self.mapping();
        let (core, begin, end) = self.parts_mut();

        SemiLocalPhmmViewMut {
            mapping,
            core,
            begin,
            end,
        }
    }
}

impl<T, const S: usize> ToOwnedData for SemiLocalPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        SemiLocalPhmm {
            mapping: self.mapping,
            core:    self.core.clone(),
            begin:   self.begin.clone(),
            end:     self.end.clone(),
        }
    }
}

impl<'a, T, const S: usize> DataView<'a> for SemiLocalPhmmView<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        SemiLocalPhmmView {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

impl<T, const S: usize> AsView for SemiLocalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        SemiLocalPhmmView {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

impl<'a, T, const S: usize> ToView<'a> for SemiLocalPhmmViewMut<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_view(self) -> Self::View<'a> {
        SemiLocalPhmmView {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

impl<T, const S: usize> ToOwnedData for SemiLocalPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        SemiLocalPhmm {
            mapping: self.mapping,
            core:    self.core.clone(),
            begin:   self.begin.clone(),
            end:     self.end.clone(),
        }
    }
}

impl<'a, T, const S: usize> DataViewMut<'a> for SemiLocalPhmmViewMut<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b, {
        SemiLocalPhmmViewMut {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

impl<T, const S: usize> AssocOwnedType for DomainPhmm<T, S> {
    type Owned = DomainPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for DomainPhmm<T, S>
where
    T: Clone + 'static,
{
    type View<'a> = DomainPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for DomainPhmm<T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = DomainPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AssocOwnedType for DomainPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type Owned = DomainPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for DomainPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type View<'a> = DomainPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for DomainPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = DomainPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AssocOwnedType for DomainPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type Owned = DomainPhmm<T, S>;
}

impl<T, const S: usize> AssocViewType for DomainPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type View<'a> = DomainPhmmView<'a, T, S>;
}

impl<T, const S: usize> AssocViewMutType for DomainPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    type ViewMut<'a> = DomainPhmmViewMut<'a, T, S>;
}

impl<T, const S: usize> AsView for DomainPhmm<T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        DomainPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   self.begin(),
            end:     self.end(),
        }
    }
}

impl<T, const S: usize> AsViewMut for DomainPhmm<T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view_mut(&mut self) -> Self::ViewMut<'_> {
        let mapping = self.mapping();
        let (core, begin, end) = self.parts_mut();

        DomainPhmmViewMut {
            mapping,
            core,
            begin,
            end,
        }
    }
}

impl<T, const S: usize> ToOwnedData for DomainPhmmView<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        DomainPhmm {
            mapping: self.mapping,
            core:    self.core.clone(),
            begin:   self.begin.clone(),
            end:     self.end.clone(),
        }
    }
}

impl<'a, T, const S: usize> DataView<'a> for DomainPhmmView<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        DomainPhmmView {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

impl<T, const S: usize> AsView for DomainPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        DomainPhmmView {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

impl<'a, T, const S: usize> ToView<'a> for DomainPhmmViewMut<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_view(self) -> Self::View<'a> {
        DomainPhmmView {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

impl<T, const S: usize> ToOwnedData for DomainPhmmViewMut<'_, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        DomainPhmm {
            mapping: self.mapping,
            core:    self.core.clone(),
            begin:   self.begin.clone(),
            end:     self.end.clone(),
        }
    }
}

impl<'a, T, const S: usize> DataViewMut<'a> for DomainPhmmViewMut<'a, T, S>
where
    T: Clone + 'static,
{
    #[inline]
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b, {
        DomainPhmmViewMut {
            mapping: self.mapping,
            core:    self.core,
            begin:   self.begin,
            end:     self.end,
        }
    }
}

/// A trait allowing a view of a global pHMM to be constructed from various
/// other pHMMs.
#[allow(dead_code)]
pub trait AsGlobalView<T, const S: usize> {
    /// Interprets the pHMM as a global pHMM via [`GlobalPhmmView`].
    #[must_use]
    fn as_global_view(&self) -> GlobalPhmmView<'_, T, S>;
}

impl<T, const S: usize> AsGlobalView<T, S> for GlobalPhmm<T, S> {
    #[inline]
    fn as_global_view(&self) -> GlobalPhmmView<'_, T, S> {
        GlobalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
        }
    }
}

impl<T, const S: usize> AsGlobalView<T, S> for DomainPhmm<T, S> {
    #[inline]
    fn as_global_view(&self) -> GlobalPhmmView<'_, T, S> {
        GlobalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
        }
    }
}

impl<T, const S: usize> AsGlobalView<T, S> for SemiLocalPhmm<T, S> {
    #[inline]
    fn as_global_view(&self) -> GlobalPhmmView<'_, T, S> {
        GlobalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
        }
    }
}

impl<T, const S: usize> AsGlobalView<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn as_global_view(&self) -> GlobalPhmmView<'_, T, S> {
        GlobalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
        }
    }
}

impl<T, const S: usize> AsGlobalView<T, S> for GlobalPhmmView<'_, T, S> {
    #[inline]
    fn as_global_view(&self) -> GlobalPhmmView<'_, T, S> {
        GlobalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
        }
    }
}

impl<T, const S: usize> AsGlobalView<T, S> for DomainPhmmView<'_, T, S> {
    #[inline]
    fn as_global_view(&self) -> GlobalPhmmView<'_, T, S> {
        GlobalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
        }
    }
}

impl<T, const S: usize> AsGlobalView<T, S> for SemiLocalPhmmView<'_, T, S> {
    #[inline]
    fn as_global_view(&self) -> GlobalPhmmView<'_, T, S> {
        GlobalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
        }
    }
}

impl<T, const S: usize> AsGlobalView<T, S> for LocalPhmmView<'_, T, S> {
    #[inline]
    fn as_global_view(&self) -> GlobalPhmmView<'_, T, S> {
        GlobalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
        }
    }
}

/// A trait allowing a view of a semilocal pHMM to be constructed from various
/// other pHMMs.
#[allow(dead_code)]
pub trait AsSemiLocalView<T, const S: usize> {
    /// Interprets the pHMM as a semilocal pHMM via [`SemiLocalPhmmView`].
    fn as_semilocal_view(&self) -> SemiLocalPhmmView<'_, T, S>;
}

impl<T, const S: usize> AsSemiLocalView<T, S> for SemiLocalPhmm<T, S> {
    #[inline]
    fn as_semilocal_view(&self) -> SemiLocalPhmmView<'_, T, S> {
        SemiLocalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   self.begin(),
            end:     self.end(),
        }
    }
}

impl<T, const S: usize> AsSemiLocalView<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn as_semilocal_view(&self) -> SemiLocalPhmmView<'_, T, S> {
        SemiLocalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   &self.begin().semilocal_params,
            end:     &self.end().semilocal_params,
        }
    }
}

impl<T, const S: usize> AsSemiLocalView<T, S> for SemiLocalPhmmView<'_, T, S> {
    #[inline]
    fn as_semilocal_view(&self) -> SemiLocalPhmmView<'_, T, S> {
        SemiLocalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   self.begin(),
            end:     self.end(),
        }
    }
}

impl<T, const S: usize> AsSemiLocalView<T, S> for LocalPhmmView<'_, T, S> {
    #[inline]
    fn as_semilocal_view(&self) -> SemiLocalPhmmView<'_, T, S> {
        SemiLocalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   &self.begin().semilocal_params,
            end:     &self.end().semilocal_params,
        }
    }
}

/// A trait allowing a view of a domain pHMM to be constructed from various
/// other pHMMs.
#[allow(dead_code)]
pub trait AsDomainView<T, const S: usize> {
    /// Interprets the pHMM as a domain pHMM via [`DomainPhmmView`].
    fn as_domain_view(&self) -> DomainPhmmView<'_, T, S>;
}

impl<T, const S: usize> AsDomainView<T, S> for DomainPhmm<T, S> {
    #[inline]
    fn as_domain_view(&self) -> DomainPhmmView<'_, T, S> {
        DomainPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   self.begin(),
            end:     self.end(),
        }
    }
}

impl<T, const S: usize> AsDomainView<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn as_domain_view(&self) -> DomainPhmmView<'_, T, S> {
        DomainPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   &self.begin().domain_params,
            end:     &self.end().domain_params,
        }
    }
}

impl<T, const S: usize> AsDomainView<T, S> for DomainPhmmView<'_, T, S> {
    #[inline]
    fn as_domain_view(&self) -> DomainPhmmView<'_, T, S> {
        DomainPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   self.begin(),
            end:     self.end(),
        }
    }
}

impl<T, const S: usize> AsDomainView<T, S> for LocalPhmmView<'_, T, S> {
    #[inline]
    fn as_domain_view(&self) -> DomainPhmmView<'_, T, S> {
        DomainPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   &self.begin().domain_params,
            end:     &self.end().domain_params,
        }
    }
}

/// A trait allowing a view of a local pHMM to be constructed from various other
/// pHMMs.
#[allow(dead_code)]
pub trait AsLocalView<T, const S: usize> {
    /// Interprets the pHMM as a local pHMM via [`LocalPhmmView`].
    fn as_local_view(&self) -> LocalPhmmView<'_, T, S>;
}

impl<T, const S: usize> AsLocalView<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn as_local_view(&self) -> LocalPhmmView<'_, T, S> {
        LocalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   self.begin(),
            end:     self.end(),
        }
    }
}

impl<T, const S: usize> AsLocalView<T, S> for LocalPhmmView<'_, T, S> {
    #[inline]
    fn as_local_view(&self) -> LocalPhmmView<'_, T, S> {
        LocalPhmmView {
            mapping: self.mapping(),
            core:    self.core(),
            begin:   self.begin(),
            end:     self.end(),
        }
    }
}

impl<T, const S: usize> PhmmIndexable for GlobalPhmmView<'_, T, S> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        self.core.num_pseudomatch()
    }
}

impl<T, const S: usize> PhmmIndexable for GlobalPhmmViewMut<'_, T, S> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        self.core.num_pseudomatch()
    }
}

impl<T, const S: usize> PhmmIndexable for DomainPhmmView<'_, T, S> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        self.core.num_pseudomatch()
    }
}

impl<T, const S: usize> PhmmIndexable for DomainPhmmViewMut<'_, T, S> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        self.core.num_pseudomatch()
    }
}

impl<T, const S: usize> PhmmIndexable for SemiLocalPhmmView<'_, T, S> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        self.core.num_pseudomatch()
    }
}

impl<T, const S: usize> PhmmIndexable for SemiLocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        self.core.num_pseudomatch()
    }
}

impl<T, const S: usize> PhmmIndexable for LocalPhmmView<'_, T, S> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        self.core.num_pseudomatch()
    }
}

impl<T, const S: usize> PhmmIndexable for LocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        self.core.num_pseudomatch()
    }
}

impl<T, const S: usize> GetCore<T, S> for GlobalPhmmView<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCore<T, S> for GlobalPhmmViewMut<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCore<T, S> for DomainPhmmView<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCore<T, S> for DomainPhmmViewMut<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCore<T, S> for SemiLocalPhmmView<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCore<T, S> for SemiLocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCore<T, S> for LocalPhmmView<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCore<T, S> for LocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCoreMut<T, S> for GlobalPhmmViewMut<'_, T, S> {
    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCoreMut<T, S> for DomainPhmmViewMut<'_, T, S> {
    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCoreMut<T, S> for SemiLocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetCoreMut<T, S> for LocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for GlobalPhmmView<'_, T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetLayer<T, S> for GlobalPhmmViewMut<'_, T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetLayer<T, S> for DomainPhmmView<'_, T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetLayer<T, S> for DomainPhmmViewMut<'_, T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetLayer<T, S> for SemiLocalPhmmView<'_, T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetLayer<T, S> for SemiLocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetLayer<T, S> for LocalPhmmView<'_, T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetLayer<T, S> for LocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for GlobalPhmmViewMut<'_, T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.core_mut().layers_mut()
    }

    #[inline]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        self.core_mut().layers_mut_vec()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for DomainPhmmViewMut<'_, T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.core_mut().layers_mut()
    }

    #[inline]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        self.core_mut().layers_mut_vec()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for SemiLocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.core_mut().layers_mut()
    }

    #[inline]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        self.core_mut().layers_mut_vec()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for LocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.core_mut().layers_mut()
    }

    #[inline]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        self.core_mut().layers_mut_vec()
    }
}

impl<T, const S: usize> GetModule for DomainPhmmView<'_, T, S> {
    type Begin = DomainModule<T, S>;
    type End = DomainModule<T, S>;

    #[inline]
    fn begin(&self) -> &Self::Begin {
        self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        self.end
    }
}

impl<T, const S: usize> GetModule for DomainPhmmViewMut<'_, T, S> {
    type Begin = DomainModule<T, S>;
    type End = DomainModule<T, S>;

    #[inline]
    fn begin(&self) -> &Self::Begin {
        self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        self.end
    }
}

impl<T, const S: usize> GetModule for SemiLocalPhmmView<'_, T, S> {
    type Begin = SemiLocalModule<T>;
    type End = SemiLocalModule<T>;

    #[inline]
    fn begin(&self) -> &Self::Begin {
        self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        self.end
    }
}

impl<T, const S: usize> GetModule for SemiLocalPhmmViewMut<'_, T, S> {
    type Begin = SemiLocalModule<T>;
    type End = SemiLocalModule<T>;

    #[inline]
    fn begin(&self) -> &Self::Begin {
        self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        self.end
    }
}

impl<T, const S: usize> GetModule for LocalPhmmView<'_, T, S> {
    type Begin = LocalModule<T, S>;
    type End = LocalModule<T, S>;

    #[inline]
    fn begin(&self) -> &Self::Begin {
        self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        self.end
    }
}

impl<T, const S: usize> GetModule for LocalPhmmViewMut<'_, T, S> {
    type Begin = LocalModule<T, S>;
    type End = LocalModule<T, S>;

    #[inline]
    fn begin(&self) -> &Self::Begin {
        self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        self.end
    }
}

impl<T, const S: usize> GetModuleMut for DomainPhmmViewMut<'_, T, S> {
    #[inline]
    fn begin_mut(&mut self) -> &mut Self::Begin {
        self.begin
    }

    #[inline]
    fn end_mut(&mut self) -> &mut Self::End {
        self.end
    }
}

impl<T, const S: usize> GetModuleMut for SemiLocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn begin_mut(&mut self) -> &mut Self::Begin {
        self.begin
    }

    #[inline]
    fn end_mut(&mut self) -> &mut Self::End {
        self.end
    }
}

impl<T, const S: usize> GetModuleMut for LocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn begin_mut(&mut self) -> &mut Self::Begin {
        self.begin
    }

    #[inline]
    fn end_mut(&mut self) -> &mut Self::End {
        self.end
    }
}

impl<T, const S: usize> GlobalPhmmView<'_, T, S> {
    /// Returns a reference to the [`ByteIndexMap`] used by the global pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GlobalPhmmViewMut<'_, T, S> {
    /// Returns a reference to the [`ByteIndexMap`] used by the global pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> LocalPhmmView<'_, T, S> {
    /// Returns a reference to the [`ByteIndexMap`] used by the local pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> LocalPhmmViewMut<'_, T, S> {
    /// Returns a reference to the [`ByteIndexMap`] used by the local pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> DomainPhmmView<'_, T, S> {
    /// Returns a reference to the [`ByteIndexMap`] used by the domain pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> DomainPhmmViewMut<'_, T, S> {
    /// Returns a reference to the [`ByteIndexMap`] used by the domain pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> SemiLocalPhmmView<'_, T, S> {
    /// Returns a reference to the [`ByteIndexMap`] used by the semilocal pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> SemiLocalPhmmViewMut<'_, T, S> {
    /// Returns a reference to the [`ByteIndexMap`] used by the semilocal pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

/// A trait providing read-only accessors to the modules at the beginning and
/// end of a pHMM view.
///
/// This is similar to [`GetModule`] but returns references with a longer
/// lifetime.
#[allow(dead_code)]
pub trait GetModuleView<'a> {
    /// The type of the module at the beginning of the pHMM.
    type Begin;
    /// The type of the module at the end of the pHMM.
    type End;

    /// Returns a reference to the module at the start of the pHMM.
    #[must_use]
    fn begin_view(&self) -> &'a Self::Begin;

    /// Returns a reference to the module at the end of the pHMM.
    #[must_use]
    fn end_view(&self) -> &'a Self::End;
}

impl<'a, T, const S: usize> GetModuleView<'a> for LocalPhmmView<'a, T, S> {
    type Begin = LocalModule<T, S>;
    type End = LocalModule<T, S>;

    #[inline]
    fn begin_view(&self) -> &'a Self::Begin {
        self.begin
    }

    #[inline]
    fn end_view(&self) -> &'a Self::End {
        self.end
    }
}

impl<'a, T, const S: usize> GetModuleView<'a> for DomainPhmmView<'a, T, S> {
    type Begin = DomainModule<T, S>;
    type End = DomainModule<T, S>;

    #[inline]
    fn begin_view(&self) -> &'a Self::Begin {
        self.begin
    }

    #[inline]
    fn end_view(&self) -> &'a Self::End {
        self.end
    }
}

impl<'a, T, const S: usize> GetModuleView<'a> for SemiLocalPhmmView<'a, T, S> {
    type Begin = SemiLocalModule<T>;
    type End = SemiLocalModule<T>;

    #[inline]
    fn begin_view(&self) -> &'a Self::Begin {
        self.begin
    }

    #[inline]
    fn end_view(&self) -> &'a Self::End {
        self.end
    }
}

/// A trait providing read-only access to the [`CorePhmm`] within a larger pHMM.
#[allow(dead_code)]
pub(crate) trait GetCoreView<'a, T, const S: usize> {
    /// Returns a reference to the [`CorePhmm`] holding the core parameters.
    #[must_use]
    fn core_view(&self) -> &'a CorePhmm<T, S>;
}

/// A trait providing read-only accessors to the layers of a pHMM view.
///
/// This is similar to [`GetLayer`] but returns references with a longer
/// lifetime.
pub trait GetLayerView<'a, T, const S: usize>: PhmmIndexable {
    /// Retrieves a slice of the layers contained within the core pHMM.
    #[must_use]
    fn layers_view(&self) -> &'a [LayerParams<T, S>];

    /// Returns the first layer, as well as all subsequent layers.
    ///
    /// This is an infallible version of `model.layers_view().split_first()`.
    #[must_use]
    fn split_first_layer_view(&self) -> (&'a LayerParams<T, S>, &'a [LayerParams<T, S>]);

    /// Returns the last layer, as well as all previous layers.
    ///
    /// This is an infallible version of `model.layers_view().split_last()`.
    #[must_use]
    fn split_last_layer_view(&self) -> (&'a LayerParams<T, S>, &'a [LayerParams<T, S>]);

    /// Returns a reference to the parameters for the layer containing the BEGIN
    /// state.
    ///
    /// This is an infallible version of `model.get_layer_view(Begin)`.
    #[inline]
    #[must_use]
    fn begin_layer_view(&self) -> &'a LayerParams<T, S> {
        &self.layers_view()[0]
    }

    /// Returns a reference to the parameters for the layer containing the first
    /// match state.
    ///
    /// This is an infallible version of `model.get_layer_view(FirstMatch)`.
    #[inline]
    #[must_use]
    fn first_match_view(&self) -> &'a LayerParams<T, S> {
        &self.layers_view()[1]
    }

    /// Returns a reference to the parameters for the layer containing the last
    /// match state.
    ///
    /// This is an infallible version of `model.get_layer_view(LastMatch)`.
    #[inline]
    #[must_use]
    fn last_match_view(&self) -> &'a LayerParams<T, S> {
        &self.layers_view()[self.layers_view().len() - 1]
    }

    /// Returns the layers, split at a particular index.
    ///
    /// If the index is out of bounds, `None` is returned.
    #[inline]
    #[must_use]
    #[allow(clippy::type_complexity)]
    fn split_layers_at_view(&self, j: impl PhmmIndex) -> Option<(&'a [LayerParams<T, S>], &'a [LayerParams<T, S>])> {
        self.layers_view().split_at_checked(self.get_dp_index(j))
    }

    /// Gets a layer from within the core pHMM.
    ///
    /// This returns `None` if the index is out of bounds or [`End`] (since
    /// there is no layer corresponding to the END state).
    ///
    /// [`End`]: crate::alignment::phmm::indexing::End
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn get_layer_view(&self, j: impl PhmmIndex) -> Option<&'a LayerParams<T, S>> {
        self.layers_view().get(self.get_dp_index(j))
    }

    /// Gets a range of layers from within the core pHMM.
    ///
    /// If any of the indices are out of bounds, this will return `None`.
    /// Particularly, if the range is end-inclusive and ends with `End` (e.g.,
    /// `..=End`), this will return `None`.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn get_layers_view(&self, range: impl PhmmIndexRange) -> Option<&'a [LayerParams<T, S>]> {
        self.layers_view().get(self.get_dp_range(range))
    }
}

impl<'a, T, const S: usize> GetCoreView<'a, T, S> for GlobalPhmmView<'a, T, S> {
    #[inline]
    fn core_view(&self) -> &'a CorePhmm<T, S> {
        self.core
    }
}

impl<'a, T, const S: usize> GetCoreView<'a, T, S> for DomainPhmmView<'a, T, S> {
    #[inline]
    fn core_view(&self) -> &'a CorePhmm<T, S> {
        self.core
    }
}

impl<'a, T, const S: usize> GetCoreView<'a, T, S> for SemiLocalPhmmView<'a, T, S> {
    #[inline]
    fn core_view(&self) -> &'a CorePhmm<T, S> {
        self.core
    }
}

impl<'a, T, const S: usize> GetCoreView<'a, T, S> for LocalPhmmView<'a, T, S> {
    #[inline]
    fn core_view(&self) -> &'a CorePhmm<T, S> {
        self.core
    }
}

impl<'a, T, const S: usize> GetLayerView<'a, T, S> for GlobalPhmmView<'a, T, S> {
    #[inline]
    fn layers_view(&self) -> &'a [LayerParams<T, S>] {
        self.core.layers()
    }

    #[inline]
    fn split_first_layer_view(&self) -> (&'a LayerParams<T, S>, &'a [LayerParams<T, S>]) {
        self.core.split_first_layer()
    }

    #[inline]
    fn split_last_layer_view(&self) -> (&'a LayerParams<T, S>, &'a [LayerParams<T, S>]) {
        self.core.split_last_layer()
    }
}

impl<'a, T, const S: usize> GetLayerView<'a, T, S> for DomainPhmmView<'a, T, S> {
    #[inline]
    fn layers_view(&self) -> &'a [LayerParams<T, S>] {
        self.core.layers()
    }

    #[inline]
    fn split_first_layer_view(&self) -> (&'a LayerParams<T, S>, &'a [LayerParams<T, S>]) {
        self.core.split_first_layer()
    }

    #[inline]
    fn split_last_layer_view(&self) -> (&'a LayerParams<T, S>, &'a [LayerParams<T, S>]) {
        self.core.split_last_layer()
    }
}

impl<'a, T, const S: usize> GetLayerView<'a, T, S> for SemiLocalPhmmView<'a, T, S> {
    #[inline]
    fn layers_view(&self) -> &'a [LayerParams<T, S>] {
        self.core.layers()
    }

    #[inline]
    fn split_first_layer_view(&self) -> (&'a LayerParams<T, S>, &'a [LayerParams<T, S>]) {
        self.core.split_first_layer()
    }

    #[inline]
    fn split_last_layer_view(&self) -> (&'a LayerParams<T, S>, &'a [LayerParams<T, S>]) {
        self.core.split_last_layer()
    }
}

impl<'a, T, const S: usize> GetLayerView<'a, T, S> for LocalPhmmView<'a, T, S> {
    #[inline]
    fn layers_view(&self) -> &'a [LayerParams<T, S>] {
        self.core.layers()
    }

    #[inline]
    fn split_first_layer_view(&self) -> (&'a LayerParams<T, S>, &'a [LayerParams<T, S>]) {
        self.core.split_first_layer()
    }

    #[inline]
    fn split_last_layer_view(&self) -> (&'a LayerParams<T, S>, &'a [LayerParams<T, S>]) {
        self.core.split_last_layer()
    }
}

impl<T, const S: usize> GetMapping<S> for GlobalPhmmView<'_, T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for GlobalPhmmViewMut<'_, T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for LocalPhmmView<'_, T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for LocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for SemiLocalPhmmView<'_, T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for SemiLocalPhmmViewMut<'_, T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for DomainPhmmView<'_, T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for DomainPhmmViewMut<'_, T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}
